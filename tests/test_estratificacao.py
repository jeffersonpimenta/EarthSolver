"""Testes do modelo direto e da estratificacao do solo."""

import numpy as np
import pytest

from earthsolver.filtros import j0, resistividade_aparente
from earthsolver.estratificacao import Estratificador
from earthsolver.solo import ModeloSolo


def _rho_a_duas_camadas_serie(rho1, rho2, h, a, n_termos=2000):
    """Forma fechada classica de Wenner para duas camadas (referencia)."""
    k = (rho2 - rho1) / (rho2 + rho1)
    n = np.arange(1, n_termos + 1)
    soma = 0.0
    for ni in n:
        soma += k ** ni * (
            1.0 / np.sqrt(1.0 + (2.0 * ni * h / a) ** 2)
            - 1.0 / np.sqrt(4.0 + (2.0 * ni * h / a) ** 2)
        )
    return rho1 * (1.0 + 4.0 * soma)


def test_j0_valores_conhecidos():
    # Zeros e valores tabelados de J0.
    assert j0(0.0) == pytest.approx(1.0, abs=1e-7)
    assert float(j0(2.4048)) == pytest.approx(0.0, abs=1e-4)  # primeiro zero
    assert float(j0(1.0)) == pytest.approx(0.7651976866, abs=1e-6)
    assert float(j0(5.0)) == pytest.approx(-0.1775967713, abs=1e-5)


def test_solo_homogeneo_devolve_rho_constante():
    # Modelo de uma camada: rho_a deve ser ~rho para qualquer espacamento.
    for a in [0.5, 2.0, 10.0, 50.0]:
        assert resistividade_aparente([300.0], [], a) == pytest.approx(300.0, rel=1e-3)


def test_duas_camadas_bate_com_serie():
    # Modelo direto (integral) vs forma fechada em series.
    rho1, rho2, h = 100.0, 500.0, 3.0
    for a in [1.0, 2.0, 5.0, 10.0, 20.0, 40.0]:
        integral = resistividade_aparente([rho1, rho2], [h], a)
        serie = _rho_a_duas_camadas_serie(rho1, rho2, h, a)
        assert integral == pytest.approx(serie, rel=2e-3)


def test_duas_camadas_baixo_para_cima():
    # Solo condutivo sobre resistivo: limites assintoticos coerentes.
    rho1, rho2, h = 80.0, 800.0, 2.0
    pequeno = resistividade_aparente([rho1, rho2], [h], 0.2)  # ve a camada 1
    grande = resistividade_aparente([rho1, rho2], [h], 200.0)  # ve a camada 2
    assert pequeno == pytest.approx(rho1, rel=0.05)
    assert grande == pytest.approx(rho2, rel=0.1)


def test_estratificar_recupera_duas_camadas():
    verdade = ModeloSolo(rho=[100.0, 400.0], espessura=[3.0])
    a = np.array([0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0])
    rho_a = resistividade_aparente(verdade.rho, verdade.espessura, a)

    estrat = Estratificador(a, resistividades=rho_a)
    modelo = estrat.estratificar(2)

    assert estrat.rms < 1.0  # ajuste praticamente perfeito (dados sinteticos)
    assert modelo.rho[0] == pytest.approx(100.0, rel=0.05)
    assert modelo.rho[1] == pytest.approx(400.0, rel=0.05)
    assert modelo.espessura[0] == pytest.approx(3.0, rel=0.1)


def test_estratificar_recupera_tres_camadas():
    verdade = ModeloSolo(rho=[300.0, 50.0, 600.0], espessura=[2.0, 5.0])
    a = np.geomspace(0.5, 60.0, 12)
    rho_a = resistividade_aparente(verdade.rho, verdade.espessura, a)

    estrat = Estratificador(a, resistividades=rho_a)
    estrat.estratificar(3)
    assert estrat.rms < 2.0


def test_auto_estratificar_escolhe_uma_camada_para_solo_uniforme():
    a = np.geomspace(0.5, 40.0, 10)
    rho_a = resistividade_aparente([250.0], [], a)
    estrat = Estratificador(a, resistividades=rho_a)
    modelo = estrat.auto_estratificar(max_camadas=3)
    assert modelo.n_camadas == 1
    assert modelo.rho[0] == pytest.approx(250.0, rel=1e-2)


def test_wenner_resistividade_a_partir_de_resistencia():
    # rho_a = 2*pi*a*R.
    estrat = Estratificador([1.0, 2.0], resistencias=[10.0, 5.0])
    assert estrat.rho_a[0] == pytest.approx(2 * np.pi * 1.0 * 10.0)
    assert estrat.rho_a[1] == pytest.approx(2 * np.pi * 2.0 * 5.0)


def test_validacoes_estratificador():
    with pytest.raises(ValueError):
        Estratificador([1.0], resistencias=[1.0])  # poucos pontos
    with pytest.raises(ValueError):
        Estratificador([2.0, 1.0], resistencias=[1.0, 1.0])  # nao crescente
    with pytest.raises(ValueError):
        Estratificador([1.0, 2.0])  # sem dados de resistencia/resistividade
