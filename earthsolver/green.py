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


def potencial_resto(rho, espessura, r_h, z, z_linha,
                    n_lambda=_N_LAMBDA, lim=_LIM):
    """Resto de camadas: G menos a parte direta+imagem (forma fechada).

    Termo limitado (finito inclusive na coincidencia fonte=observacao), obtido
    pela integral de Hankel do problema de contorno por lambda. Zero para solo
    uniforme. Aceita arrays (broadcast); retorna a mesma forma achatada nao e
    necessaria pois quem chama (montagem da matriz) ja trabalha vetorizado.
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

    interfaces = np.cumsum(espessura)            # profundidades d_1..d_{N-1}
    K = 2 * N - 1                                # incognitas: A_j,B_j (j<N-1) + A_{N-1}

    def iA(j):
        return 2 * j if j < N - 1 else 2 * N - 2

    # --- grade de lambda (resto decai ~ e^{-lambda * escala de profundidade}) ---
    ell = max(min(espessura.min(), float((zz + zl).min())), 1e-3)
    lam = np.linspace(0.0, lim / ell, n_lambda)
    lam[0] = lam[1] * 1e-3                        # evita matriz singular em lambda=0

    eh = np.exp(-lam[:, None] * espessura[None, :])   # (L, N-1): e^{-lambda h_j}

    # --- matriz do sistema (independe de z, z'): montar e inverter uma vez ---
    M = np.zeros((lam.size, K, K))
    M[:, 0, iA(0)] = -1.0                         # superficie: -A_0 + B_0 e^{-lh_0} = ...
    M[:, 0, 1] = eh[:, 0]
    for i in range(N - 1):                        # interface i (camadas i, i+1)
        rp, rc = 1 + 2 * i, 2 + 2 * i            # linhas: continuidade psi / corrente
        M[:, rp, iA(i)] += eh[:, i]
        M[:, rp, 2 * i + 1] += 1.0
        M[:, rp, iA(i + 1)] += -1.0
        if i + 1 < N - 1:
            M[:, rp, 2 * (i + 1) + 1] += -eh[:, i + 1]
        M[:, rc, iA(i)] += -eh[:, i] / rho[i]
        M[:, rc, 2 * i + 1] += 1.0 / rho[i]
        M[:, rc, iA(i + 1)] += 1.0 / rho[i + 1]
        if i + 1 < N - 1:
            M[:, rc, 2 * (i + 1) + 1] += -eh[:, i + 1] / rho[i + 1]
    Minv = np.linalg.inv(M)                       # (L, K, K)

    out = np.empty(zz.size, dtype=float)
    for e in range(zz.size):
        rh_e, z_e, zl_e = rh[e], zz[e], zl[e]
        m = _camada(zl_e, interfaces)            # camada da fonte
        n = _camada(z_e, interfaces)             # camada da observacao
        prefm = rho[m] / (4.0 * np.pi)

        # termo primario (fonte) movido para o RHS: M x = -primario
        c = np.zeros((lam.size, K))
        if m == 0:
            c[:, 0] = -prefm * np.exp(-lam * zl_e)
        for i in range(N - 1):
            sgn = (1.0 if m == i else 0.0) - (1.0 if m == i + 1 else 0.0)
            if sgn == 0.0:
                continue
            di = interfaces[i]
            decai = np.exp(-lam * abs(di - zl_e))
            # continuidade de psi usa o primario direto (prefm = rho_m/4pi);
            # continuidade de (1/rho) dpsi/dz cancela o rho_m -> fator 1/(4pi).
            c[:, 1 + 2 * i] += -sgn * prefm * decai
            c[:, 2 + 2 * i] += sgn * np.sign(di - zl_e) * decai / (4.0 * np.pi)

        x = np.einsum("lij,lj->li", Minv, c)     # coeficientes A_j, B_j por lambda

        # potencial homogeneo na camada de observacao
        topo_n = 0.0 if n == 0 else interfaces[n - 1]
        psi_h = x[:, iA(n)] * np.exp(-lam * (z_e - topo_n))
        if n < N - 1:
            psi_h = psi_h + x[:, 2 * n + 1] * np.exp(-lam * (interfaces[n] - z_e))
        prim = prefm * np.exp(-lam * abs(z_e - zl_e)) if n == m else 0.0

        fechada = prefm * (np.exp(-lam * abs(z_e - zl_e))
                           + np.exp(-lam * (z_e + zl_e)))
        bracket = psi_h + prim - fechada
        out[e] = _trapz(bracket * j0(lam * rh_e), lam)

    return float(out[0]) if forma == () else out.reshape(forma)
