"""Funcao de Green do solo estratificado em N camadas.

Potencial de uma fonte pontual de corrente enterrada num solo horizontalmente
estratificado, com a interface ar/solo (imagem de superficie). Base do solver de
segmentacao (numerico.py).

Formulacao: para cada numero de onda lambda, o potencial transformado em Hankel
satisfaz um problema de contorno 1-D em z. Em cada camada
    psi_j(z) = A_j e^{-lambda (z - topo_j)} + B_j e^{-lambda (base_j - z)}
(expoentes <= 0 dentro da camada -> sem overflow). As condicoes de contorno sao:
continuidade de psi e de (1/rho) dpsi/dz nas interfaces, dpsi/dz = 0 na superficie
(ar isolante -> imagem) e decaimento no semi-espaco inferior. Montamos o sistema
linear (2N-1 incognitas) por lambda e resolvemos. O potencial fisico e a integral
de Hankel
    G(r_h, z, z') = integral_0^inf psi_obs(lambda; z, z') J0(lambda r_h) dlambda.

Para condicionar a integral, subtraimos as formas fechadas direta + imagem
(transformadas exatas 1/r e 1/r_img) e integramos so o resto limitado, reusando
j0() de filtros.py (mesmo espirito do truque de convergencia de filtros.py).
O resultado e exato para qualquer par de camadas (a subtracao apenas melhora o
condicionamento; ver deducao no design).

Caso uniforme (N=1): forma fechada rho/(4pi) (1/r + 1/r_img), sem integral.
"""

import numpy as np

from .filtros import j0

# np.trapz foi renomeado para np.trapezoid no numpy 2.0.
_trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))

# Profundidade adimensional ate a qual o resto e desprezivel (e^-12 ~ 6e-6 por
# reflexao) e numero de pontos da quadratura de Hankel.
_LIM = 12.0
_N_LAMBDA = 4001


def _camada(prof, interfaces):
    """Indice (0-based) da camada que contem a profundidade `prof`.

    interfaces = profundidades acumuladas das N-1 interfaces (cumsum(espessura)).
    """
    return int(np.searchsorted(interfaces, prof, side="right"))


def potencial_camadas(rho, espessura, r_h, z, z_linha,
                      n_lambda=_N_LAMBDA, lim=_LIM):
    """Potencial de uma fonte pontual de corrente unitaria no solo N-camadas.

    rho: resistividades das camadas (rho_1..rho_N), Ohm.m.
    espessura: espessuras das N-1 camadas superiores (a ultima e infinita), m.
    r_h: distancia horizontal fonte->observacao (m).
    z: profundidade da observacao (m, >= 0).
    z_linha: profundidade da fonte (m, >= 0).
    Aceita r_h, z, z_linha escalares ou arrays (broadcast). Retorna G (V por
    corrente unitaria, i.e. resistencia mutua de ponto), mesma forma da entrada.

    G = forma fechada (direta + imagem de superficie) + resto de camadas, este
    ultimo dado por potencial_resto (limitado). Reusa a parte fechada exata.
    """
    rho = np.asarray(rho, dtype=float)
    espessura = np.asarray(espessura, dtype=float)
    N = rho.size

    forma = np.broadcast_shapes(np.shape(r_h), np.shape(z), np.shape(z_linha))
    rh = np.broadcast_to(np.asarray(r_h, dtype=float), forma).ravel()
    zz = np.broadcast_to(np.asarray(z, dtype=float), forma).ravel()
    zl = np.broadcast_to(np.asarray(z_linha, dtype=float), forma).ravel()

    interfaces = np.cumsum(espessura) if N > 1 else np.array([])
    m = (np.searchsorted(interfaces, zl, side="right") if N > 1
         else np.zeros(zz.size, dtype=int))
    prefm = rho[m] / (4.0 * np.pi)
    r = np.sqrt(rh ** 2 + (zz - zl) ** 2)
    rimg = np.sqrt(rh ** 2 + (zz + zl) ** 2)
    out = prefm * (1.0 / r + 1.0 / rimg)
    if N > 1:
        out = out + potencial_resto(rho, espessura, rh, zz, zl,
                                    n_lambda=n_lambda, lim=lim)
    return float(out[0]) if forma == () else out.reshape(forma)


def preparar_resto(rho, espessura, ell, n_lambda=_N_LAMBDA, lim=_LIM):
    """Precomputa o sistema do problema de contorno por lambda (independe dos
    pontos): grade `lam`, `interfaces` e a inversa `Minv` do sistema (L, K, K).

    Reusar entre muitas avaliacoes do resto evita reconstruir/inverter o sistema
    a cada chamada -- era um custo repetido O(M*ng) vezes na montagem da matriz.
    `ell` = escala de profundidade (m) que fixa o alcance de lambda.
    """
    rho = np.asarray(rho, dtype=float)
    espessura = np.asarray(espessura, dtype=float)
    N = rho.size
    interfaces = np.cumsum(espessura)            # profundidades d_1..d_{N-1}
    K = 2 * N - 1                                # incognitas: A_j,B_j (j<N-1) + A_{N-1}

    def iA(j):
        return 2 * j if j < N - 1 else 2 * N - 2

    ell = max(float(ell), 1e-3)
    lam = np.linspace(0.0, lim / ell, n_lambda)
    lam[0] = lam[1] * 1e-3                        # evita matriz singular em lambda=0
    eh = np.exp(-lam[:, None] * espessura[None, :])   # (L, N-1): e^{-lambda h_j}

    Mat = np.zeros((lam.size, K, K))
    Mat[:, 0, iA(0)] = -1.0                       # superficie: -A_0 + B_0 e^{-lh_0} = ...
    Mat[:, 0, 1] = eh[:, 0]
    for i in range(N - 1):                        # interface i (camadas i, i+1)
        rp, rc = 1 + 2 * i, 2 + 2 * i            # linhas: continuidade psi / corrente
        Mat[:, rp, iA(i)] += eh[:, i]
        Mat[:, rp, 2 * i + 1] += 1.0
        Mat[:, rp, iA(i + 1)] += -1.0
        if i + 1 < N - 1:
            Mat[:, rp, 2 * (i + 1) + 1] += -eh[:, i + 1]
        Mat[:, rc, iA(i)] += -eh[:, i] / rho[i]
        Mat[:, rc, 2 * i + 1] += 1.0 / rho[i]
        Mat[:, rc, iA(i + 1)] += 1.0 / rho[i + 1]
        if i + 1 < N - 1:
            Mat[:, rc, 2 * (i + 1) + 1] += -eh[:, i + 1] / rho[i + 1]
    Minv = np.linalg.inv(Mat)                     # (L, K, K)
    return {"rho": rho, "espessura": espessura, "interfaces": interfaces,
            "N": N, "K": K, "lam": lam, "Minv": Minv}


def _topo(interfaces, n):
    """Profundidade do topo da camada n (0 para a 1a camada). n array (P,)."""
    idx = np.clip(n - 1, 0, max(interfaces.size - 1, 0))
    return np.where(n > 0, interfaces[idx], 0.0)


def brackets_resto(ctx, z, zl):
    """`bracket(lambda)` (psi_h + primario - forma fechada) para pares de
    profundidades (observacao z, fonte zl), arrays (P,). Depende so das
    profundidades, nao de r_h. Retorna (P, L). Vetorizado sobre P.
    """
    rho = ctx["rho"]; interfaces = ctx["interfaces"]
    N = ctx["N"]; K = ctx["K"]; lam = ctx["lam"]; Minv = ctx["Minv"]
    L = lam.size
    z = np.asarray(z, dtype=float).ravel()
    zl = np.asarray(zl, dtype=float).ravel()
    P = z.size

    msrc = np.searchsorted(interfaces, zl, side="right")    # camada da fonte (P,)
    nobs = np.searchsorted(interfaces, z, side="right")     # camada da obs   (P,)
    prefm = rho[msrc] / (4.0 * np.pi)                        # (P,)

    # primario movido ao RHS: Mat x = -primario  -> c (P, L, K)
    c = np.zeros((P, L, K))
    is0 = msrc == 0
    if is0.any():
        c[is0, :, 0] = -prefm[is0, None] * np.exp(-lam[None, :] * zl[is0, None])
    for i in range(N - 1):
        sgn = (msrc == i).astype(float) - (msrc == i + 1).astype(float)   # (P,)
        nz = sgn != 0.0
        if not nz.any():
            continue
        di = interfaces[i]
        decai = np.exp(-lam[None, :] * np.abs(di - zl[nz, None]))          # (Pnz, L)
        c[nz, :, 1 + 2 * i] += (-sgn[nz, None] * prefm[nz, None]) * decai
        c[nz, :, 2 + 2 * i] += (sgn[nz, None] * np.sign(di - zl[nz, None])
                                * decai) / (4.0 * np.pi)

    x = np.einsum("lij,plj->pli", Minv, c)                  # coef. por lambda (P, L, K)

    p_ar = np.arange(P)
    iAn = np.where(nobs < N - 1, 2 * nobs, 2 * N - 2)
    psi_h = x[p_ar, :, iAn] * np.exp(-lam[None, :]
                                     * (z[:, None] - _topo(interfaces, nobs)[:, None]))
    mask = nobs < N - 1
    if mask.any():
        idx2 = np.clip(2 * nobs + 1, 0, K - 1)
        base_n = interfaces[np.clip(nobs, 0, max(N - 2, 0))]
        extra = x[p_ar, :, idx2] * np.exp(-lam[None, :] * (base_n[:, None] - z[:, None]))
        psi_h = psi_h + np.where(mask[:, None], extra, 0.0)

    dz = np.abs(z - zl)[:, None]
    prim = np.where((nobs == msrc)[:, None],
                    prefm[:, None] * np.exp(-lam[None, :] * dz), 0.0)
    fechada = prefm[:, None] * (np.exp(-lam[None, :] * dz)
                                + np.exp(-lam[None, :] * (z + zl)[:, None]))
    return psi_h + prim - fechada


def integrar_hankel(ctx, brackets, idx, rh, chunk=8192):
    """Integral de Hankel vetorizada por blocos: para cada ponto,
    trapz(brackets[idx] * J0(lam * rh), lam). brackets (U, L); idx (P,); rh (P,).
    `chunk` limita o array transitorio (chunk, L). Retorna (P,).
    """
    lam = ctx["lam"]
    rh = np.asarray(rh, dtype=float).ravel()
    idx = np.asarray(idx).ravel()
    out = np.empty(rh.size, dtype=float)
    for s in range(0, rh.size, chunk):
        e = min(s + chunk, rh.size)
        out[s:e] = _trapz(brackets[idx[s:e]] * j0(lam[None, :] * rh[s:e, None]),
                          lam, axis=1)
    return out


def grade_rh(rh_max, ell, n_lin=400, n_log=200):
    """Grade de distancias horizontais para tabular o resto g(r_h).

    O resto e suave em r_h e varia na escala ~`ell` (decai p/ r_h >> ell), entao
    densifica perto de 0 (linear ate ~6*ell) e usa log ate `rh_max`. Para
    geometria puramente vertical (rh_max=0) devolve so [0.0] (g avaliado em r_h=0,
    exato). Reutilizada por _resto_camadas e _superficie_resto.
    """
    if rh_max <= 0.0:
        return np.array([0.0])
    rt = min(max(6.0 * float(ell), 1e-3), rh_max)
    g = np.linspace(0.0, rt, n_lin)
    if rh_max > rt:
        g = np.concatenate([g, np.geomspace(rt, rh_max, n_log)])
    return np.unique(g)


def tabela_hankel(ctx, brk, npair, rh_grid):
    """Tabela g(r_h) por par de profundidades: integral de Hankel de cada bracket
    nos nos `rh_grid`. brk (npair, L) -> (npair, G). Custo O(npair*G*L), feito UMA
    vez; depois a montagem so interpola (sem reintegrar por distancia)."""
    G = rh_grid.size
    idx = np.repeat(np.arange(npair), G)
    rh = np.tile(rh_grid, npair)
    return integrar_hankel(ctx, brk, idx, rh).reshape(npair, G)


def interp_tabela(table, rh_grid, pair, rh):
    """Interpolacao linear vetorizada da tabela do resto.

    table (npair, G) amostrada em `rh_grid` (G,). `pair` (indice de par de
    profundidade) e `rh` tem a mesma forma -> devolve g dessa forma. Usa indexacao
    avancada nos dois eixos (sem laco Python sobre pares)."""
    G = rh_grid.size
    pair = np.asarray(pair)
    if G == 1:                                   # so r_h=0 (geometria vertical)
        return table[pair, 0]
    k = np.clip(np.searchsorted(rh_grid, rh) - 1, 0, G - 2)
    x0 = rh_grid[k]
    t = (rh - x0) / (rh_grid[k + 1] - x0)
    return table[pair, k] * (1.0 - t) + table[pair, k + 1] * t


def avaliador_resto(ctx, brk, npair, rh_grid, n_direct_total):
    """Escolhe a estrategia mais barata p/ avaliar o resto g(par, r_h).

    - tabela + interpolacao: ~`npair*G` integrais de Hankel, independente de M;
      vence quando ha muitos pares (matriz M^2, raster grande).
    - integral direta por chamada: ~`n_direct_total` integrais; vence em problemas
      pequenos (poucos pares), onde tabular sairia mais caro que so integrar.
    Devolve um callable evalg(pair, rh) -> g, com a forma de `rh`."""
    if npair * rh_grid.size <= n_direct_total:
        table = tabela_hankel(ctx, brk, npair, rh_grid)

        def evalg(pair, rh):
            return interp_tabela(table, rh_grid, pair, rh)
    else:
        def evalg(pair, rh):
            rh = np.asarray(rh)
            return integrar_hankel(ctx, brk, np.asarray(pair).ravel(),
                                  rh.ravel()).reshape(rh.shape)
    return evalg


def potencial_resto(rho, espessura, r_h, z, z_linha,
                    n_lambda=_N_LAMBDA, lim=_LIM, ctx=None):
    """Resto de camadas: G menos a parte direta+imagem (forma fechada).

    Termo limitado (finito inclusive na coincidencia fonte=observacao), obtido
    pela integral de Hankel do problema de contorno por lambda. Zero para solo
    uniforme. Aceita arrays (broadcast). Vetorizado: deduplica pares de
    profundidade (z, z') -> brackets e integra J0 em blocos. `ctx` (de
    preparar_resto) pode ser passado para reusar o sistema entre chamadas.
    """
    rho = np.asarray(rho, dtype=float)
    espessura = np.asarray(espessura, dtype=float)
    N = rho.size

    forma = np.broadcast_shapes(np.shape(r_h), np.shape(z), np.shape(z_linha))
    rh = np.broadcast_to(np.asarray(r_h, dtype=float), forma).ravel()
    zz = np.broadcast_to(np.asarray(z, dtype=float), forma).ravel()
    zl = np.broadcast_to(np.asarray(z_linha, dtype=float), forma).ravel()
    if N == 1:
        out = np.zeros(zz.size)
        return out[0] if forma == () else out.reshape(forma)

    if ctx is None:
        ell = min(espessura.min(), float((zz + zl).min()))
        ctx = preparar_resto(rho, espessura, ell, n_lambda=n_lambda, lim=lim)

    pares = np.stack([zz, zl], axis=1)
    uniq, inv = np.unique(pares, axis=0, return_inverse=True)
    bracket = brackets_resto(ctx, uniq[:, 0], uniq[:, 1])     # (U, L)
    out = integrar_hankel(ctx, bracket, inv.ravel(), rh)
    return float(out[0]) if forma == () else out.reshape(forma)
