"""Solver numerico de segmentacao de condutores (geometria).

Condutor: um segmento reto de eletrodo (p1->p2) com raio. Eletrodo: a colecao de
condutores, com segmentacao automatica (segmentar) e um construtor de malha
retangular regular (malha_retangular) para reproduzir o exemplo do IEEE Std 80.
A montagem da matriz e a solucao ficam em EstudoNumerico (proxima fase).
"""

import math
from dataclasses import dataclass

import numpy as np

from .green import potencial_resto
from .seguranca import tensoes_toleraveis
from .solo import ModeloSolo


@dataclass
class Condutor:
    """Segmento reto de eletrodo: p1, p2 = (x, y, z) em metros, raio em metros.

    Convencao: z = profundidade >= 0 (superficie em z=0, condutor enterrado).
    """

    p1: tuple
    p2: tuple
    raio: float

    def __post_init__(self):
        self.p1 = tuple(float(v) for v in self.p1)
        self.p2 = tuple(float(v) for v in self.p2)
        self.raio = float(self.raio)
        if self.raio <= 0:
            raise ValueError("raio deve ser > 0")
        if self.p1[2] < 0 or self.p2[2] < 0:
            raise ValueError("z deve ser >= 0 (condutor enterrado)")
        if self.comprimento <= 0:
            raise ValueError("comprimento do condutor deve ser > 0")

    @property
    def comprimento(self) -> float:
        return math.dist(self.p1, self.p2)


@dataclass
class Segmentos:
    """Segmentos discretizados (arrays vetorizados de M segmentos).

    mid (M,3): pontos medios; dir (M,3): versores; meia (M,): meio-comprimento;
    raio (M,): raios; camada (M,): indice 0-based da camada de solo do segmento.
    """

    mid: np.ndarray
    dir: np.ndarray
    meia: np.ndarray
    raio: np.ndarray
    camada: np.ndarray

    @property
    def n(self) -> int:
        return len(self.meia)


class Eletrodo:
    """Colecao de condutores de um sistema de aterramento."""

    def __init__(self, condutores):
        self.condutores = list(condutores)
        if not self.condutores:
            raise ValueError("eletrodo precisa de ao menos um condutor")

    def segmentar(self, comp_alvo, interfaces=()) -> Segmentos:
        """Subdivide cada condutor em segmentos de comprimento <= comp_alvo.

        Quebra os condutores nas profundidades `interfaces` (limites de camada)
        para que cada segmento fique inteiramente em uma unica camada de solo.
        """
        if comp_alvo <= 0:
            raise ValueError("comp_alvo deve ser > 0")
        interfaces = np.asarray(interfaces, dtype=float)
        mids, dirs, meias, raios, camadas = [], [], [], [], []
        for c in self.condutores:
            p1, p2 = np.array(c.p1), np.array(c.p2)
            L = c.comprimento
            u = (p2 - p1) / L
            # pontos de quebra (parametro t em [0,1]) nas interfaces cruzadas
            ts = {0.0, 1.0}
            dz = p2[2] - p1[2]
            if abs(dz) > 1e-12 and interfaces.size:
                for d_int in interfaces:
                    t = (d_int - p1[2]) / dz
                    if 1e-9 < t < 1.0 - 1e-9:
                        ts.add(float(t))
            ts = sorted(ts)
            for ta, tb in zip(ts[:-1], ts[1:]):
                sub = L * (tb - ta)
                nseg = max(1, int(math.ceil(sub / comp_alvo - 1e-9)))
                for k in range(nseg):
                    a = p1 + (p2 - p1) * (ta + (tb - ta) * k / nseg)
                    b = p1 + (p2 - p1) * (ta + (tb - ta) * (k + 1) / nseg)
                    mid = 0.5 * (a + b)
                    mids.append(mid)
                    dirs.append(u)
                    meias.append(0.5 * float(np.linalg.norm(b - a)))
                    raios.append(c.raio)
                    camadas.append(int(np.searchsorted(interfaces, mid[2], side="right")))
        return Segmentos(np.array(mids), np.array(dirs), np.array(meias),
                         np.array(raios), np.array(camadas, dtype=int))

    @classmethod
    def malha_retangular(cls, comprimento_x, comprimento_y, espac_x, espac_y,
                         prof_h, d, n_hastes=0, comp_haste=0.0):
        """Gera o reticulado retangular regular (na profundidade prof_h) + hastes.

        Linhas paralelas a x (e a y) espacadas de espac_y (e espac_x). Hastes
        verticais distribuidas igualmente nos nos do perimetro, descendo de
        prof_h ate prof_h+comp_haste. Raio = d/2.
        """
        raio = d / 2.0
        nx = int(round(comprimento_x / espac_x))
        ny = int(round(comprimento_y / espac_y))
        xs = np.linspace(0.0, comprimento_x, nx + 1)
        ys = np.linspace(0.0, comprimento_y, ny + 1)

        condutores = []
        for y in ys:                                   # linhas paralelas a x
            condutores.append(Condutor((0.0, y, prof_h),
                                       (comprimento_x, y, prof_h), raio))
        for x in xs:                                   # linhas paralelas a y
            condutores.append(Condutor((x, 0.0, prof_h),
                                       (x, comprimento_y, prof_h), raio))

        if n_hastes > 0:
            if comp_haste <= 0:
                raise ValueError("comp_haste deve ser > 0 quando ha hastes")
            loop = cls._nos_perimetro(xs, ys)
            if n_hastes > len(loop):
                raise ValueError(
                    f"n_hastes ({n_hastes}) excede os nos do perimetro "
                    f"({len(loop)}); hastes coincidiriam")
            for i in range(n_hastes):
                x, y = loop[int(round(i * len(loop) / n_hastes)) % len(loop)]
                condutores.append(Condutor((x, y, prof_h),
                                           (x, y, prof_h + comp_haste), raio))
        return cls(condutores)

    @staticmethod
    def _nos_perimetro(xs, ys):
        """Nos da grade no perimetro, em ordem de laco (sem repetir cantos)."""
        loop = [(x, ys[0]) for x in xs]
        loop += [(xs[-1], y) for y in ys[1:]]
        loop += [(x, ys[-1]) for x in xs[-2::-1]]
        loop += [(xs[0], y) for y in ys[-2:0:-1]]
        return loop


def _gauss(n):
    """Nos (em [-1,1]) e pesos de Gauss-Legendre normalizados (soma dos pesos = 1)."""
    xi, w = np.polynomial.legendre.leggauss(n)
    return xi, w / 2.0


def _phi_seg(P, A, B, rho, raio):
    """Potencial em P (Nf,3) de cada segmento-fonte A->B (M,3) com corrente
    unitaria uniforme, kernel rho/(4 pi r). Retorna (Nf, M).

    Forma fechada da integral de linha (arcsinh):
        integral_0^L dl/|P - (A + l u)| = asinh((L - s0)/p) + asinh(s0/p),
    com p = distancia perpendicular. A distancia perpendicular e regularizada
    pelo raio do condutor (p_eff = sqrt(perp^2 + a^2)), o que dá o valor correto
    tambem para segmentos colineares (vizinhos na mesma reta -> perp = 0) e e a
    mesma regularizacao que origina o termo proprio ln(2L/a).
    """
    AB = B - A
    L = np.linalg.norm(AB, axis=1)
    u = AB / L[:, None]
    PA = P[:, None, :] - A[None, :, :]
    d1sq = np.einsum("nmk,nmk->nm", PA, PA)
    s0 = np.einsum("nmk,mk->nm", PA, u)
    pe = np.sqrt(np.maximum(d1sq - s0 ** 2, 0.0) + raio[None, :] ** 2)
    val = np.arcsinh((L[None, :] - s0) / pe) + np.arcsinh(s0 / pe)
    return rho[None, :] / (4.0 * np.pi * L[None, :]) * val


class EstudoNumerico:
    """Estudo de aterramento pelo metodo numerico de segmentacao de condutores.

    Discretiza o eletrodo em segmentos, monta a matriz de resistencias R (potencial
    medio em cada segmento por corrente drenada pelos demais) usando a funcao de
    Green do solo N-camadas (green.py), resolve o sistema equipotencial para a
    distribuicao de corrente e devolve Rg, GPR e (proxima fase) tensoes de toque/passo.

    modelo_solo: ModeloSolo (N camadas). eletrodo: Eletrodo. Ig: corrente injetada (A).
    t: duracao da falta (s). peso: 50 ou 70 (kg). comp_alvo: comprimento alvo de
    segmento (m). rho_s/h_s: camada superficial (brita) para tensoes toleraveis.
    n_gauss: ordem da quadratura externa por segmento.
    """

    def __init__(self, modelo_solo, eletrodo, Ig, t, peso=70, comp_alvo=2.0,
                 rho_s=None, h_s=0.1, n_gauss=4, passo_raster=2.0, margem_raster=5.0):
        if not isinstance(modelo_solo, ModeloSolo):
            raise ValueError("modelo_solo deve ser uma instancia de ModeloSolo")
        if not isinstance(eletrodo, Eletrodo):
            raise ValueError("eletrodo deve ser uma instancia de Eletrodo")
        if Ig is None or Ig <= 0:
            raise ValueError("Ig deve ser > 0")
        if t is None or t <= 0:
            raise ValueError("t deve ser > 0")
        if peso not in (50, 70):
            raise ValueError("peso deve ser 50 ou 70 (kg)")
        if rho_s is not None and rho_s <= 0:
            raise ValueError("rho_s deve ser > 0")
        if h_s <= 0:
            raise ValueError("h_s deve ser > 0")
        if n_gauss < 1:
            raise ValueError("n_gauss deve ser >= 1")

        self.solo = modelo_solo
        self.eletrodo = eletrodo
        self.Ig = float(Ig)
        self.t = float(t)
        self.peso = peso
        self.n_gauss = n_gauss
        self.passo_raster = float(passo_raster)
        self.margem_raster = float(margem_raster)

        interfaces = np.cumsum(modelo_solo.espessura)
        self.segs = eletrodo.segmentar(comp_alvo, interfaces)
        self.rho_seg = np.array([modelo_solo.rho[c] for c in self.segs.camada])
        self.rho_eq = modelo_solo.uniforme_equivalente()
        self.rho_s = float(rho_s) if rho_s is not None else self.rho_eq
        self.h_s = float(h_s)

        self.R = None
        self.I = None
        self.V = None
        self.Em = None
        self.Es = None
        self.raster = None
        self.raster_toque = None
        self.raster_passo = None
        self.resultado = None

    # --------------------------------------------------------- matriz R
    def montar_resistencias(self):
        """Monta a matriz de resistencias R (M x M), simetrica positiva-definida."""
        s = self.segs
        M = s.n
        ng = self.n_gauss
        meia, dir_, mid = s.meia, s.dir, s.mid
        raio = s.raio
        L = 2.0 * meia
        A = mid - meia[:, None] * dir_
        B = mid + meia[:, None] * dir_
        esp = np.array([0.0, 0.0, 1.0])                 # espelho z -> -z
        A_img = A - 2.0 * A[:, 2][:, None] * esp
        B_img = B - 2.0 * B[:, 2][:, None] * esp
        self._A, self._B, self._A_img, self._B_img = A, B, A_img, B_img

        xi, wg = _gauss(ng)
        # pontos de Gauss de cada segmento (campo): (M, ng, 3)
        Pg = mid[:, None, :] + xi[None, :, None] * meia[:, None, None] * dir_[:, None, :]
        P = Pg.reshape(M * ng, 3)

        phi_dir = _phi_seg(P, A, B, self.rho_seg, raio)        # (Nf, M)
        phi_img = _phi_seg(P, A_img, B_img, self.rho_seg, raio)
        phi = (phi_dir + phi_img).reshape(M, ng, M)
        R = np.einsum("a,iaj->ij", wg, phi)

        # diagonal: auto-potencial direto (forma fechada) + imagem (quadratura)
        self_dir = self.rho_seg / (2.0 * np.pi * L) * (np.log(2.0 * L / raio) - 1.0)
        self_img = np.einsum("a,iai->i", wg, phi_img.reshape(M, ng, M))
        np.fill_diagonal(R, self_dir + self_img)

        if self.solo.n_camadas > 1:
            R = R + self._resto_camadas(Pg, wg)

        self.R = 0.5 * (R + R.T)                          # simetriza
        return self.R

    def _resto_camadas(self, Pg, wg):
        """Termo de correcao de camadas (G - direta-imagem), por quadratura dupla.

        So entra quando o solo tem mais de uma camada; para solo uniforme e zero.
        """
        s = self.segs
        M = s.n
        ng = self.n_gauss
        Pf = Pg.reshape(M * ng, 3)
        Rrem = np.zeros((M, M))
        for j in range(M):                               # fonte: segmento j
            for b in range(ng):
                xs, ys, zs = Pg[j, b]
                rh = np.hypot(Pf[:, 0] - xs, Pf[:, 1] - ys)
                grem = potencial_resto(self.solo.rho, self.solo.espessura,
                                       rh, Pf[:, 2], zs).reshape(M, ng)
                Rrem[:, j] += wg[b] * (grem @ wg)
        return Rrem

    # --------------------------------------------------------- solucao
    def resolver_rg(self) -> float:
        """Resolve so a resistencia de malha Rg (sem o campo de superficie).

        Monta R (se preciso), resolve o sistema equipotencial e fixa V (GPR) e I.
        E o caminho rapido usado por resolver() e pela curva de convergencia.
        """
        if self.R is None:
            self.montar_resistencias()
        y = np.linalg.solve(self.R, np.ones(self.segs.n))
        Rg = 1.0 / float(y.sum())
        self.V = self.Ig * Rg
        self.I = self.V * y
        return Rg

    def resolver(self) -> dict:
        """Resolve o sistema equipotencial e o estudo completo de seguranca.

        Devolve um dicionario com chaves compativeis com o metodo IEEE 80
        (Rg, GPR, Em, Es, E_toque, E_passo, *_ok, aprovado) mais extras
        (n_segmentos, V, rho_eq). O raster do potencial fica em self.raster.
        """
        Rg = self.resolver_rg()

        self._calcular_superficie()
        E_toque, E_passo, Cs = tensoes_toleraveis(self.rho_eq, self.rho_s,
                                                  self.h_s, self.t, self.peso)
        self.resultado = {
            "Rg": Rg,
            "GPR": self.V,
            "V": self.V,
            "Em": self.Em,
            "Es": self.Es,
            "E_toque": E_toque,
            "E_passo": E_passo,
            "Cs": Cs,
            "toque_ok": bool(self.Em <= E_toque),
            "passo_ok": bool(self.Es <= E_passo),
            "aprovado": bool(self.Em <= E_toque and self.Es <= E_passo),
            "n_segmentos": self.segs.n,
            "rho_eq": self.rho_eq,
        }
        return self.resultado

    # --------------------------------------------------------- potencial de superficie
    def potencial_superficie(self, xy):
        """Potencial na superficie (z=0) nos pontos xy (lista/array de (x,y)).

        Phi(p) = soma_j I_j * G(p, segmento_j). Requer resolver() previo (I).
        """
        if self.I is None:
            raise ValueError("rode resolver() antes de consultar o potencial")
        xy = np.atleast_2d(np.asarray(xy, dtype=float))
        P = np.column_stack([xy[:, 0], xy[:, 1], np.zeros(len(xy))])
        raio = self.segs.raio
        phi = (_phi_seg(P, self._A, self._B, self.rho_seg, raio)
               + _phi_seg(P, self._A_img, self._B_img, self.rho_seg, raio))
        Phi = phi @ self.I
        if self.solo.n_camadas > 1:
            Phi = Phi + self._superficie_resto(P)
        return Phi

    def _superficie_resto(self, P):
        """Correcao de camadas no potencial de superficie (solo N-camadas)."""
        s = self.segs
        ng = self.n_gauss
        xi, wg = _gauss(ng)
        mid, meia, dir_ = s.mid, s.meia, s.dir
        Pg = mid[:, None, :] + xi[None, :, None] * meia[:, None, None] * dir_[:, None, :]
        extra = np.zeros(len(P))
        for j in range(s.n):
            for b in range(ng):
                xs, ys, zs = Pg[j, b]
                rh = np.hypot(P[:, 0] - xs, P[:, 1] - ys)
                grem = potencial_resto(self.solo.rho, self.solo.espessura,
                                       rh, P[:, 2], zs)
                extra += wg[b] * self.I[j] * grem
        return extra

    def _calcular_superficie(self):
        """Monta o raster de potencial e extrai tensoes de toque (Em) e passo (Es)."""
        pts = np.vstack([self._A[:, :2], self._B[:, :2]])
        xmin, ymin = pts.min(axis=0)
        xmax, ymax = pts.max(axis=0)
        passo, margem = self.passo_raster, self.margem_raster
        xs = np.arange(xmin - margem, xmax + margem + 1e-9, passo)
        ys = np.arange(ymin - margem, ymax + margem + 1e-9, passo)
        X, Y = np.meshgrid(xs, ys)
        base = np.column_stack([X.ravel(), Y.ravel()])
        Phi = self.potencial_superficie(base)
        self.raster = (X, Y, Phi.reshape(X.shape))

        # campo de toque: GPR - potencial de superficie (todo o raster).
        self.raster_toque = (X, Y, self.V - Phi.reshape(X.shape))

        # toque: pior caso = GPR - menor potencial dentro da projecao da malha.
        dentro = (base[:, 0] >= xmin) & (base[:, 0] <= xmax) \
            & (base[:, 1] >= ymin) & (base[:, 1] <= ymax)
        if not dentro.any():                       # eletrodo degenerado (ex.: haste)
            dentro = np.ones(len(base), dtype=bool)
        self.Em = self.V - float(Phi[dentro].min())

        # campo de passo: maior diferenca de potencial entre pontos a 1 m
        # (em x e em y), por ponto. Es = pior caso de todo o campo.
        dphi_x = self.potencial_superficie(base + [1.0, 0.0]) - Phi
        dphi_y = self.potencial_superficie(base + [0.0, 1.0]) - Phi
        passo_pt = np.maximum(np.abs(dphi_x), np.abs(dphi_y))
        self.raster_passo = (X, Y, passo_pt.reshape(X.shape))
        self.Es = float(passo_pt.max())

    def dados_corrente(self):
        """Geometria e corrente drenada por segmento: (A, B, I).

        A, B (M,3): extremos de cada segmento; I (M,): corrente drenada para o
        solo por segmento (picos nos cantos da malha). Requer resolver() previo.
        """
        if self.I is None:
            raise ValueError("rode resolver() antes de consultar a corrente")
        return self._A, self._B, self.I

    # --------------------------------------------------------- relatorio
    def imprimir_resultado(self, cd: int = 2) -> None:
        if self.resultado is None:
            self.resolver()
        r = self.resultado
        print("Estudo de Aterramento (metodo numerico de segmentacao):")
        print("-" * 56)
        print(f"  Segmentos                 : {r['n_segmentos']}")
        print(f"  Resistencia de malha  Rg  : {r['Rg']:.{cd}f} Ohm")
        print(f"  Elevacao de potencial GPR : {r['GPR']:.{cd}f} V")
        print("-" * 56)
        print(f"  Tensao de toque  Em       : {r['Em']:.{cd}f} V "
              f"(limite {r['E_toque']:.{cd}f} V) -> "
              f"{'OK' if r['toque_ok'] else 'NAO ATENDE'}")
        print(f"  Tensao de passo  Es       : {r['Es']:.{cd}f} V "
              f"(limite {r['E_passo']:.{cd}f} V) -> "
              f"{'OK' if r['passo_ok'] else 'NAO ATENDE'}")
        print("-" * 56)
        print(f"  Veredito: {'APROVADO' if r['aprovado'] else 'REPROVADO'} "
              f"(peso {self.peso} kg, t = {self.t} s)")
        print("-" * 56)

    def exportar(self, arquivo: str = "aterramento_numerico.json") -> None:
        import json
        if self.resultado is None:
            self.resolver()
        dados = {
            "entrada": {
                "Ig": self.Ig, "t": self.t, "peso": self.peso,
                "metodo": "numerico", "rho_s": self.rho_s, "h_s": self.h_s,
                "n_condutores": len(self.eletrodo.condutores),
                "solo": self.solo.to_dict(),
            },
            "resultado": self.resultado,
        }
        with open(arquivo, "w") as f:
            json.dump(dados, f, indent=4)
        print(f"Resultado exportado para {arquivo}.")

    def exportar_raster(self, arquivo: str = "potencial.json") -> None:
        """Exporta o mapa de potencial de superficie (x, y, Phi) para plotagem."""
        import json
        if self.raster is None:
            self.resolver()
        X, Y, Phi = self.raster
        dados = {"x": X.tolist(), "y": Y.tolist(), "phi": Phi.tolist(),
                 "GPR": self.V}
        with open(arquivo, "w") as f:
            json.dump(dados, f)
        print(f"Raster de potencial exportado para {arquivo}.")


def estudo_convergencia(modelo_solo, eletrodo, Ig, t, comp_alvos, **kwargs):
    """Curva de convergencia: Rg e GPR por numero de segmentos.

    Resolve o eletrodo para cada comprimento-alvo em `comp_alvos` (caminho rapido
    resolver_rg, sem o campo de superficie) e devolve os resultados ordenados por
    numero crescente de segmentos. `kwargs` extras vao para EstudoNumerico (peso,
    rho_s, h_s, n_gauss, ...). Devolve {n_segmentos, Rg, GPR} (arrays).
    """
    kwargs.pop("comp_alvo", None)
    ns, rgs, gprs = [], [], []
    for ca in comp_alvos:
        est = EstudoNumerico(modelo_solo, eletrodo, Ig, t, comp_alvo=ca, **kwargs)
        rg = est.resolver_rg()
        ns.append(est.segs.n)
        rgs.append(rg)
        gprs.append(est.V)
    ordem = np.argsort(ns)
    return {"n_segmentos": np.asarray(ns)[ordem],
            "Rg": np.asarray(rgs)[ordem],
            "GPR": np.asarray(gprs)[ordem]}

