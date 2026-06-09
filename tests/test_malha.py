"""Testes da simulacao de malha (IEEE Std 80).

O cenario de referencia e a malha quadrada do exemplo classico do IEEE Std 80
(solo uniforme rho=400 Ohm.m, 70 m x 70 m, espacamento 7 m, profundidade 0.5 m),
para a qual a resistencia de Sverak e ~2.78 Ohm.
"""

import math

import pytest

from earthsolver.malha import EstudoAterramento, Malha
from earthsolver.solo import ModeloSolo


def _malha_exemplo(**kw):
    base = dict(
        area=70.0 * 70.0,
        Lc=11 * 70.0 * 2,        # 11 condutores em cada direcao, 70 m
        comprimento_x=70.0,
        comprimento_y=70.0,
        espac_D=7.0,
        prof_h=0.5,
        d=0.01,
    )
    base.update(kw)
    return Malha(**base)


def test_resistencia_sverak_exemplo_ieee():
    solo = ModeloSolo(rho=[400.0], espessura=[])
    estudo = EstudoAterramento(solo, _malha_exemplo(), Ig=1908.0, t=0.5, peso=70)
    r = estudo.resolver()
    assert r["Rg"] == pytest.approx(2.78, abs=0.05)


def test_fator_n_malha_quadrada():
    solo = ModeloSolo(rho=[400.0], espessura=[])
    estudo = EstudoAterramento(solo, _malha_exemplo(), Ig=1908.0, t=0.5)
    # Malha quadrada: n = na = 2*Lc/perimetro = 2*1540/280 = 11.
    assert estudo._fator_n() == pytest.approx(11.0, rel=1e-6)


def test_gpr_e_relacoes_basicas():
    solo = ModeloSolo(rho=[400.0], espessura=[])
    estudo = EstudoAterramento(solo, _malha_exemplo(), Ig=1908.0, t=0.5)
    r = estudo.resolver()
    assert r["GPR"] == pytest.approx(1908.0 * r["Rg"])
    assert r["Em"] > 0 and r["Es"] > 0
    assert r["Em"] < r["GPR"]  # tensao de toque e fracao do GPR


def test_tensoes_toleraveis_formula():
    # Cenario com brita: rho=400, rho_s=2500, h_s=0.102, t=0.5, 70 kg.
    solo = ModeloSolo(rho=[400.0], espessura=[])
    malha = _malha_exemplo(rho_s=2500.0, h_s=0.102)
    estudo = EstudoAterramento(solo, malha, Ig=1908.0, t=0.5, peso=70)
    E_toque, E_passo, Cs = estudo.tensoes_toleraveis()

    Cs_ref = 1.0 - 0.09 * (1.0 - 400.0 / 2500.0) / (2.0 * 0.102 + 0.09)
    assert Cs == pytest.approx(Cs_ref, rel=1e-9)
    c = 0.157
    assert E_toque == pytest.approx((1000 + 1.5 * Cs_ref * 2500) * c / math.sqrt(0.5), rel=1e-9)
    assert E_passo == pytest.approx((1000 + 6.0 * Cs_ref * 2500) * c / math.sqrt(0.5), rel=1e-9)


def test_cs_unitario_sem_brita():
    solo = ModeloSolo(rho=[400.0], espessura=[])
    estudo = EstudoAterramento(solo, _malha_exemplo(), Ig=1000.0, t=0.5)
    _, _, Cs = estudo.tensoes_toleraveis()
    assert Cs == pytest.approx(1.0, abs=1e-9)  # rho_s = rho


def test_peso_50_mais_restritivo_que_70():
    solo = ModeloSolo(rho=[400.0], espessura=[])
    malha = _malha_exemplo(rho_s=2500.0, h_s=0.102)
    e50 = EstudoAterramento(solo, malha, Ig=1908.0, t=0.5, peso=50)
    e70 = EstudoAterramento(solo, malha, Ig=1908.0, t=0.5, peso=70)
    t50, p50, _ = e50.tensoes_toleraveis()
    t70, p70, _ = e70.tensoes_toleraveis()
    assert t50 < t70 and p50 < p70


def test_hastes_reduzem_resistencia():
    solo = ModeloSolo(rho=[400.0], espessura=[])
    sem = EstudoAterramento(solo, _malha_exemplo(), Ig=1908.0, t=0.5).resolver()
    com = EstudoAterramento(
        solo, _malha_exemplo(n_hastes=20, comp_haste=7.5), Ig=1908.0, t=0.5
    ).resolver()
    assert com["Rg"] < sem["Rg"]


def test_validacoes_malha_e_estudo():
    with pytest.raises(ValueError):
        Malha(area=-1, Lc=1, comprimento_x=1, comprimento_y=1,
              espac_D=1, prof_h=1, d=0.01)
    solo = ModeloSolo(rho=[400.0], espessura=[])
    with pytest.raises(ValueError):
        EstudoAterramento(solo, _malha_exemplo(), Ig=0, t=0.5)
    with pytest.raises(ValueError):
        EstudoAterramento(solo, _malha_exemplo(), Ig=100, t=0.5, peso=60)
    with pytest.raises(ValueError, match="EstudoNumerico"):
        EstudoAterramento(solo, _malha_exemplo(), Ig=100, t=0.5, metodo="numerico")
