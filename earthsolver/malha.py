"""Simulacao de malha de aterramento pelo metodo simplificado do IEEE Std 80.

Malha descreve a geometria do reticulado; EstudoAterramento calcula a resistencia
de aterramento (Rg, Sverak), a elevacao de potencial (GPR) e as tensoes de malha
(toque) e de passo, comparando-as com os limites toleraveis de seguranca
(IEEE Std 80 / NBR 15751).
"""

import json
import math
from dataclasses import dataclass

from .seguranca import C_PESO, fator_cs
from .seguranca import tensoes_toleraveis as _tensoes_toleraveis
from .solo import ModeloSolo

# Profundidade de referencia para o fator Kh (IEEE 80), em metros.
_H0 = 1.0


@dataclass
class Malha:
    """Geometria de uma malha de aterramento retangular.

    area: area coberta pela malha (m^2).
    Lc: comprimento total de condutores horizontais (m).
    comprimento_x, comprimento_y: dimensoes do retangulo (m).
    espac_D: espacamento entre condutores paralelos (m).
    prof_h: profundidade de enterramento da malha (m).
    d: diametro do condutor (m).
    n_hastes: numero de hastes de aterramento.
    comp_haste: comprimento de cada haste (m).
    rho_s: resistividade da camada superficial (brita), Ohm.m. Se None, usa o
        rho do solo (sem camada de brita, Cs = 1).
    h_s: espessura da camada superficial (m).
    """

    area: float
    Lc: float
    comprimento_x: float
    comprimento_y: float
    espac_D: float
    prof_h: float
    d: float
    n_hastes: int = 0
    comp_haste: float = 0.0
    rho_s: float = None
    h_s: float = 0.1

    def __post_init__(self):
        positivos = {
            "area": self.area, "Lc": self.Lc,
            "comprimento_x": self.comprimento_x,
            "comprimento_y": self.comprimento_y,
            "espac_D": self.espac_D, "prof_h": self.prof_h, "d": self.d,
        }
        for nome, val in positivos.items():
            if val is None or val <= 0:
                raise ValueError(f"{nome} deve ser > 0 (recebido {val!r})")
        if self.n_hastes < 0:
            raise ValueError("n_hastes deve ser >= 0")
        if self.n_hastes > 0 and self.comp_haste <= 0:
            raise ValueError("comp_haste deve ser > 0 quando ha hastes")
        if self.rho_s is not None and self.rho_s <= 0:
            raise ValueError("rho_s deve ser > 0")
        if self.h_s <= 0:
            raise ValueError("h_s deve ser > 0")

    @property
    def Lr(self) -> float:
        """Comprimento total de hastes (m)."""
        return self.n_hastes * self.comp_haste

    @property
    def L_T(self) -> float:
        """Comprimento total enterrado: condutores + hastes (m)."""
        return self.Lc + self.Lr

    @property
    def perimetro(self) -> float:
        return 2.0 * (self.comprimento_x + self.comprimento_y)


class EstudoAterramento:
    """Estudo de aterramento de uma malha pelo metodo do IEEE Std 80.

    Parametros:
      modelo_solo: ModeloSolo (reduzido a uniforme equivalente para as formulas).
      malha: objeto Malha com a geometria.
      Ig: corrente de malha (A) que efetivamente escoa para o solo.
      t: duracao do choque / da falta (s).
      peso: peso corporeo de referencia, 50 ou 70 (kg).
      metodo: "ieee80" (atual). Reservado para "numerico" no futuro.
    """

    _C_PESO = C_PESO

    def __init__(self, modelo_solo, malha, Ig, t, peso=70, metodo="ieee80"):
        if not isinstance(modelo_solo, ModeloSolo):
            raise ValueError("modelo_solo deve ser uma instancia de ModeloSolo")
        if not isinstance(malha, Malha):
            raise ValueError("malha deve ser uma instancia de Malha")
        if Ig is None or Ig <= 0:
            raise ValueError("Ig deve ser > 0")
        if t is None or t <= 0:
            raise ValueError("t deve ser > 0")
        if peso not in self._C_PESO:
            raise ValueError("peso deve ser 50 ou 70 (kg)")
        if metodo == "numerico":
            raise ValueError(
                "o metodo 'numerico' (segmentacao de condutores) usa a classe "
                "EstudoNumerico, que recebe geometria explicita de condutores "
                "(earthsolver.numerico), nao a Malha agregada do IEEE 80"
            )
        if metodo != "ieee80":
            raise ValueError(f"metodo {metodo!r} desconhecido (use 'ieee80')")

        self.modelo_solo = modelo_solo
        self.malha = malha
        self.Ig = float(Ig)
        self.t = float(t)
        self.peso = peso
        self.metodo = metodo

        self.rho = modelo_solo.uniforme_equivalente()
        self.rho_s = malha.rho_s if malha.rho_s is not None else self.rho
        self.resultado = None

    # ----------------------------------------------------------- fatores IEEE80
    def _fator_n(self) -> float:
        """Fator geometrico n = na*nb*nc*nd (IEEE 80)."""
        m = self.malha
        na = 2.0 * m.Lc / m.perimetro
        # nb=1 para malha quadrada; geral usa a razao perimetro/area.
        nb = math.sqrt(m.perimetro / (4.0 * math.sqrt(m.area)))
        nc = 1.0  # retangular
        nd = 1.0  # retangular
        return na * nb * nc * nd

    def _Kii(self, n: float) -> float:
        """Fator de correcao de condutores internos (1 se ha hastes)."""
        if self.malha.n_hastes > 0:
            return 1.0
        return 1.0 / (2.0 * n) ** (2.0 / n)

    def _Km(self, n: float) -> float:
        m = self.malha
        D, h, d = m.espac_D, m.prof_h, m.d
        Kii = self._Kii(n)
        Kh = math.sqrt(1.0 + h / _H0)
        termo1 = math.log(
            D ** 2 / (16.0 * h * d)
            + (D + 2.0 * h) ** 2 / (8.0 * D * d)
            - h / (4.0 * d)
        )
        termo2 = (Kii / Kh) * math.log(8.0 / (math.pi * (2.0 * n - 1.0)))
        return (termo1 + termo2) / (2.0 * math.pi)

    def _Ks(self, n: float) -> float:
        m = self.malha
        D, h = m.espac_D, m.prof_h
        return (1.0 / math.pi) * (
            1.0 / (2.0 * h)
            + 1.0 / (D + h)
            + (1.0 / D) * (1.0 - 0.5 ** (n - 2.0))
        )

    def _resistencia_sverak(self) -> float:
        """Resistencia de aterramento Rg pela formula de Sverak (IEEE 80)."""
        m = self.malha
        A, h, L_T = m.area, m.prof_h, m.L_T
        return self.rho * (
            1.0 / L_T
            + 1.0 / math.sqrt(20.0 * A)
            * (1.0 + 1.0 / (1.0 + h * math.sqrt(20.0 / A)))
        )

    def _comprimentos_efetivos(self, diagonal: float):
        """Comprimentos efetivos L_M (malha/toque) e L_S (passo)."""
        m = self.malha
        if m.n_hastes > 0:
            fator = 1.55 + 1.22 * (m.comp_haste / diagonal)
            L_M = m.Lc + fator * m.Lr
        else:
            L_M = m.Lc + m.Lr
        L_S = 0.75 * m.Lc + 0.85 * m.Lr
        return L_M, L_S

    # ----------------------------------------------------------- seguranca
    def _Cs(self) -> float:
        """Fator de reducao da camada superficial (IEEE 80)."""
        return fator_cs(self.rho, self.rho_s, self.malha.h_s)

    def tensoes_toleraveis(self):
        """Tensoes de toque e passo toleraveis (V) para o peso de referencia."""
        return _tensoes_toleraveis(self.rho, self.rho_s, self.malha.h_s,
                                   self.t, self.peso)

    # ----------------------------------------------------------- solucao
    def resolver(self) -> dict:
        """Executa o estudo completo e devolve um dicionario de resultados."""
        m = self.malha
        diagonal = math.sqrt(m.comprimento_x ** 2 + m.comprimento_y ** 2)
        n = self._fator_n()
        Ki = 0.644 + 0.148 * n
        Km = self._Km(n)
        Ks = self._Ks(n)
        L_M, L_S = self._comprimentos_efetivos(diagonal)

        Rg = self._resistencia_sverak()
        GPR = self.Ig * Rg
        Em = self.rho * self.Ig * Km * Ki / L_M     # tensao de malha (toque)
        Es = self.rho * self.Ig * Ks * Ki / L_S     # tensao de passo

        E_toque, E_passo, Cs = self.tensoes_toleraveis()

        self.resultado = {
            "rho": self.rho,
            "rho_s": self.rho_s,
            "n": n,
            "Ki": Ki,
            "Km": Km,
            "Ks": Ks,
            "Kh": math.sqrt(1.0 + m.prof_h / _H0),
            "Kii": self._Kii(n),
            "Cs": Cs,
            "L_M": L_M,
            "L_S": L_S,
            "Rg": Rg,
            "GPR": GPR,
            "Em": Em,
            "Es": Es,
            "E_toque": E_toque,
            "E_passo": E_passo,
            "toque_ok": Em <= E_toque,
            "passo_ok": Es <= E_passo,
            "aprovado": (Em <= E_toque) and (Es <= E_passo),
        }
        return self.resultado

    # ----------------------------------------------------------- relatorio
    def imprimir_resultado(self, cd: int = 2) -> None:
        if self.resultado is None:
            self.resolver()
        r = self.resultado
        print("Estudo de Aterramento (IEEE Std 80):")
        print("-" * 48)
        print(f"  Resistividade equivalente : {r['rho']:.{cd}f} Ohm.m")
        print(f"  Resistencia de malha  Rg  : {r['Rg']:.{cd}f} Ohm")
        print(f"  Elevacao de potencial GPR : {r['GPR']:.{cd}f} V")
        print("-" * 48)
        print(f"  Tensao de toque  Em       : {r['Em']:.{cd}f} V "
              f"(limite {r['E_toque']:.{cd}f} V) -> "
              f"{'OK' if r['toque_ok'] else 'NAO ATENDE'}")
        print(f"  Tensao de passo  Es       : {r['Es']:.{cd}f} V "
              f"(limite {r['E_passo']:.{cd}f} V) -> "
              f"{'OK' if r['passo_ok'] else 'NAO ATENDE'}")
        print("-" * 48)
        print(f"  Veredito: {'APROVADO' if r['aprovado'] else 'REPROVADO'} "
              f"(peso {self.peso} kg, t = {self.t} s)")
        print("-" * 48)

    def exportar(self, arquivo: str = "aterramento.json") -> None:
        if self.resultado is None:
            self.resolver()
        dados = {
            "entrada": {
                "Ig": self.Ig, "t": self.t, "peso": self.peso,
                "metodo": self.metodo, "malha": vars(self.malha),
                "solo": self.modelo_solo.to_dict(),
            },
            "resultado": self.resultado,
        }
        with open(arquivo, "w") as f:
            json.dump(dados, f, indent=4)
        print(f"Resultado exportado para {arquivo}.")
