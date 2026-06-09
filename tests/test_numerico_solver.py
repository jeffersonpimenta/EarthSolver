"""Testes do nucleo do solver numerico: matriz de resistencias, solucao, Rg.

Oraculos:
  - haste vertical unica vs formula analitica R = rho/(2 pi L)(ln(8L/d) - 1);
  - matriz R simetrica e definida positiva;
  - exemplo IEEE Std 80 (70x70, solo uniforme 400) -> Rg ~ 2.78 Ohm (Sverak);
  - hastes reduzem Rg; convergencia de Rg ao refinar a segmentacao.
"""

import math

import numpy as np
import pytest

from earthsolver.numerico import Condutor, Eletrodo, EstudoNumerico
from earthsolver.solo import ModeloSolo


def _malha_ieee(**kw):
    return Eletrodo.malha_retangular(70.0, 70.0, 7.0, 7.0, 0.5, 0.01, **kw)


def test_haste_vertical_vs_analitico():
    solo = ModeloSolo([100.0], [])
    haste = Eletrodo([Condutor((0, 0, 0.0), (0, 0, 3.0), 0.01)])
    r = EstudoNumerico(solo, haste, Ig=10.0, t=0.5, comp_alvo=0.25).resolver()
    L, d, rho = 3.0, 0.02, 100.0
    Ra = rho / (2.0 * math.pi * L) * (math.log(8.0 * L / d) - 1.0)
    assert r["Rg"] == pytest.approx(Ra, rel=0.05)


def test_matriz_simetrica_definida_positiva():
    solo = ModeloSolo([400.0], [])
    est = EstudoNumerico(solo, _malha_ieee(), Ig=1000.0, t=0.5, comp_alvo=7.0)
    est.resolver()
    R = est.R
    assert np.allclose(R, R.T, rtol=1e-9, atol=1e-9)
    assert np.all(np.linalg.eigvalsh(R) > 0.0)          # definida positiva


def test_ieee80_rg_proximo_sverak():
    solo = ModeloSolo([400.0], [])
    r = EstudoNumerico(solo, _malha_ieee(), Ig=1908.0, t=0.5, comp_alvo=3.5).resolver()
    assert r["Rg"] == pytest.approx(2.78, rel=0.08)


def test_gpr_consistente():
    solo = ModeloSolo([400.0], [])
    r = EstudoNumerico(solo, _malha_ieee(), Ig=1908.0, t=0.5, comp_alvo=7.0).resolver()
    assert r["GPR"] == pytest.approx(1908.0 * r["Rg"])
    assert r["n_segmentos"] > 0


def test_hastes_reduzem_rg():
    solo = ModeloSolo([400.0], [])
    sem = EstudoNumerico(solo, _malha_ieee(), Ig=1908.0, t=0.5, comp_alvo=7.0).resolver()
    com = EstudoNumerico(solo, _malha_ieee(n_hastes=20, comp_haste=7.5),
                         Ig=1908.0, t=0.5, comp_alvo=7.0).resolver()
    assert com["Rg"] < sem["Rg"]


def test_convergencia_rg_refinando():
    solo = ModeloSolo([400.0], [])
    grosso = EstudoNumerico(solo, _malha_ieee(), Ig=1908.0, t=0.5, comp_alvo=7.0).resolver()
    fino = EstudoNumerico(solo, _malha_ieee(), Ig=1908.0, t=0.5, comp_alvo=3.5).resolver()
    assert fino["Rg"] == pytest.approx(grosso["Rg"], rel=0.05)
