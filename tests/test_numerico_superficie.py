"""Testes do campo de potencial de superficie, toque/passo e raster (numerico)."""

import numpy as np
import pytest

from earthsolver.numerico import Eletrodo, EstudoNumerico
from earthsolver.solo import ModeloSolo


def _est(rho_s=None, h_s=0.1, n_hastes=0, comp_haste=0.0):
    solo = ModeloSolo([400.0], [])
    el = Eletrodo.malha_retangular(70.0, 70.0, 7.0, 7.0, 0.5, 0.01,
                                   n_hastes=n_hastes, comp_haste=comp_haste)
    return EstudoNumerico(solo, el, Ig=1908.0, t=0.5, peso=70,
                          comp_alvo=7.0, rho_s=rho_s, h_s=h_s)


def test_resultado_chaves_ieee_compativeis():
    r = _est(rho_s=2500.0, h_s=0.102).resolver()
    for k in ("Rg", "GPR", "Em", "Es", "E_toque", "E_passo",
              "toque_ok", "passo_ok", "aprovado", "n_segmentos"):
        assert k in r


def test_toque_passo_positivos_menores_que_gpr():
    r = _est(rho_s=2500.0, h_s=0.102).resolver()
    assert r["Em"] > 0 and r["Es"] > 0
    assert r["Em"] < r["GPR"]
    assert isinstance(r["aprovado"], bool)


def test_potencial_superficie_decai_longe():
    est = _est()
    est.resolver()
    centro = est.potencial_superficie([[35.0, 35.0]])[0]      # sobre a malha
    longe = est.potencial_superficie([[2000.0, 2000.0]])[0]   # muito distante
    assert centro > longe
    assert centro < est.V                                     # < GPR
    assert longe < 0.05 * est.V                               # decai longe


def test_raster_finito_e_consistente():
    est = _est()
    est.resolver()
    X, Y, Phi = est.raster
    assert Phi.shape == X.shape == Y.shape
    assert np.all(np.isfinite(Phi))
    assert Phi.max() <= est.V + 1e-6                          # nada acima do GPR


def test_brita_eleva_tensoes_toleraveis():
    sem = _est().resolver()                                    # rho_s = rho -> Cs = 1
    com = _est(rho_s=2500.0, h_s=0.102).resolver()
    assert com["E_toque"] > sem["E_toque"]
    assert com["E_passo"] > sem["E_passo"]
