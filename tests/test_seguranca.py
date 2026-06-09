"""Testes das tensoes toleraveis compartilhadas (seguranca.py)."""

import math

import pytest

from earthsolver.seguranca import C_PESO, fator_cs, tensoes_toleraveis


def test_cs_formula():
    Cs = fator_cs(rho=400.0, rho_s=2500.0, h_s=0.102)
    esperado = 1.0 - 0.09 * (1.0 - 400.0 / 2500.0) / (2.0 * 0.102 + 0.09)
    assert Cs == pytest.approx(esperado, rel=1e-12)


def test_cs_unitario_sem_brita():
    assert fator_cs(rho=400.0, rho_s=400.0, h_s=0.1) == pytest.approx(1.0, abs=1e-12)


def test_tensoes_toleraveis_peso70():
    E_t, E_s, Cs = tensoes_toleraveis(rho=400.0, rho_s=2500.0, h_s=0.102,
                                      t=0.5, peso=70)
    Cs_ref = 1.0 - 0.09 * (1.0 - 400.0 / 2500.0) / (2.0 * 0.102 + 0.09)
    c = C_PESO[70]
    assert Cs == pytest.approx(Cs_ref, rel=1e-12)
    assert E_t == pytest.approx((1000 + 1.5 * Cs_ref * 2500) * c / math.sqrt(0.5), rel=1e-12)
    assert E_s == pytest.approx((1000 + 6.0 * Cs_ref * 2500) * c / math.sqrt(0.5), rel=1e-12)


def test_peso_50_mais_restritivo():
    t50, p50, _ = tensoes_toleraveis(400.0, 2500.0, 0.102, 0.5, peso=50)
    t70, p70, _ = tensoes_toleraveis(400.0, 2500.0, 0.102, 0.5, peso=70)
    assert t50 < t70 and p50 < p70


def test_peso_invalido():
    with pytest.raises(ValueError):
        tensoes_toleraveis(400.0, 2500.0, 0.102, 0.5, peso=60)
