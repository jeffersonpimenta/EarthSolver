"""Testes da funcao de Green do solo N-camadas (green.py).

Oraculos independentes que fixam a formula antes do solver:
  1. Limite uniforme (1 camada) -> fonte + imagem de superficie (forma fechada).
  2. Finitude em r_h -> 0 (termo proprio).
  3. Reciprocidade G(z,z') == G(z',z).
  4. Duas camadas vs serie de imagens de Tagg/Sunde (coef. k).
"""

import math

import pytest

from earthsolver.green import potencial_camadas


def _g_uniforme(rho, r_h, z, zl):
    """Forma fechada do semi-espaco uniforme: rho/(4pi) (1/r + 1/r_imagem)."""
    r = math.sqrt(r_h ** 2 + (z - zl) ** 2)
    rimg = math.sqrt(r_h ** 2 + (z + zl) ** 2)
    return rho / (4.0 * math.pi) * (1.0 / r + 1.0 / rimg)


def _g_duas_camadas_serie(rho1, rho2, h, r_h, z, zl, n_termos=400):
    """Oraculo de duas camadas (fonte e observacao na camada 1) pela serie de
    imagens: G = rho1/(4pi) * sum_{n=-inf}^{inf} k^|n| [1/r1n + 1/r2n],
    com k=(rho2-rho1)/(rho2+rho1),
        r1n = sqrt(r_h^2 + (z - zl - 2 n h)^2),
        r2n = sqrt(r_h^2 + (z + zl - 2 n h)^2).
    Reduz ao uniforme quando k=0 (so n=0 sobrevive) e satisfaz Neumann em z=0.
    """
    k = (rho2 - rho1) / (rho2 + rho1)
    total = 0.0
    for n in range(-n_termos, n_termos + 1):
        kn = k ** abs(n)
        if kn == 0.0:
            continue
        r1n = math.sqrt(r_h ** 2 + (z - zl - 2.0 * n * h) ** 2)
        r2n = math.sqrt(r_h ** 2 + (z + zl - 2.0 * n * h) ** 2)
        total += kn * (1.0 / r1n + 1.0 / r2n)
    return rho1 / (4.0 * math.pi) * total


def test_uniforme_reduz_a_fonte_mais_imagem():
    rho = 100.0
    G = potencial_camadas([rho], [], r_h=3.0, z=0.5, z_linha=0.5)
    assert G == pytest.approx(_g_uniforme(rho, 3.0, 0.5, 0.5), rel=1e-6)


def test_uniforme_rh_zero_finito():
    rho = 100.0
    G = potencial_camadas([rho], [], r_h=0.0, z=0.6, z_linha=0.5)
    assert math.isfinite(G)
    assert G == pytest.approx(_g_uniforme(rho, 0.0, 0.6, 0.5), rel=1e-6)


def test_reciprocidade_duas_camadas():
    rho, esp = [100.0, 300.0], [2.0]
    G1 = potencial_camadas(rho, esp, r_h=4.0, z=0.5, z_linha=1.2)
    G2 = potencial_camadas(rho, esp, r_h=4.0, z=1.2, z_linha=0.5)
    assert G1 == pytest.approx(G2, rel=1e-6)


@pytest.mark.parametrize("rho2", [300.0, 30.0])
def test_duas_camadas_vs_serie_imagens(rho2):
    rho1, h = 100.0, 2.0
    z, zl, r_h = 0.5, 0.8, 5.0          # ambos na camada 1 (z, zl < h)
    G = potencial_camadas([rho1, rho2], [h], r_h=r_h, z=z, z_linha=zl)
    esperado = _g_duas_camadas_serie(rho1, rho2, h, r_h, z, zl)
    assert G == pytest.approx(esperado, rel=1e-4)


def test_reciprocidade_fonte_em_camada_profunda():
    # observacao na camada 1 (z=0.5), fonte na camada 2 (zl=3.0 > h=2): cruza camada.
    rho, esp = [100.0, 300.0], [2.0]
    G1 = potencial_camadas(rho, esp, r_h=3.0, z=0.5, z_linha=3.0)
    G2 = potencial_camadas(rho, esp, r_h=3.0, z=3.0, z_linha=0.5)
    assert G1 == pytest.approx(G2, rel=1e-5)


def test_tres_camadas_iguais_reduz_a_uniforme():
    # rho identico nas 3 camadas -> deve coincidir com o semi-espaco uniforme.
    G3 = potencial_camadas([200.0, 200.0, 200.0], [1.5, 2.5],
                           r_h=4.0, z=0.7, z_linha=1.3)
    assert G3 == pytest.approx(_g_uniforme(200.0, 4.0, 0.7, 1.3), rel=1e-5)
