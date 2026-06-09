"""Testes do solver numerico em solo estratificado (caminho do resto de camadas).

Geometria pequena (haste atravessando a interface) p/ exercitar a correcao de
camadas com custo baixo.
"""

import numpy as np
import pytest

from earthsolver.numerico import Condutor, Eletrodo, EstudoNumerico
from earthsolver.solo import ModeloSolo


def _haste(L=3.0):
    return Eletrodo([Condutor((0.0, 0.0, 0.0), (0.0, 0.0, L), 0.01)])


def test_duas_camadas_iguais_reduz_uniforme():
    # rho identico nas duas camadas -> a correcao de camadas deve ser ~0.
    uni = EstudoNumerico(ModeloSolo([200.0], []), _haste(),
                         Ig=10.0, t=0.5, comp_alvo=0.5).resolver()
    dup = EstudoNumerico(ModeloSolo([200.0, 200.0], [1.5]), _haste(),
                         Ig=10.0, t=0.5, comp_alvo=0.5).resolver()
    assert dup["Rg"] == pytest.approx(uni["Rg"], rel=2e-3)


def test_fundo_resistivo_aumenta_rg():
    uni = EstudoNumerico(ModeloSolo([100.0], []), _haste(),
                         Ig=10.0, t=0.5, comp_alvo=0.5).resolver()
    dois = EstudoNumerico(ModeloSolo([100.0, 1000.0], [1.5]), _haste(),
                          Ig=10.0, t=0.5, comp_alvo=0.5).resolver()
    assert dois["Rg"] > uni["Rg"]                        # camada de fundo resistiva


def test_matriz_camadas_simetrica_definida_positiva():
    est = EstudoNumerico(ModeloSolo([100.0, 400.0], [1.5]), _haste(),
                         Ig=10.0, t=0.5, comp_alvo=0.5)
    est.resolver()
    R = est.R
    assert np.allclose(R, R.T, atol=1e-9)
    assert np.all(np.linalg.eigvalsh(R) > 0.0)
