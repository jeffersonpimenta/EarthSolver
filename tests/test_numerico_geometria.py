"""Testes de geometria do solver numerico (Condutor, Eletrodo, segmentar)."""


import numpy as np
import pytest

from earthsolver.numerico import Condutor, Eletrodo


def test_condutor_comprimento():
    c = Condutor((0.0, 0.0, 0.5), (3.0, 4.0, 0.5), 0.005)
    assert c.comprimento == pytest.approx(5.0)


def test_condutor_validacao():
    with pytest.raises(ValueError):
        Condutor((0, 0, 0.5), (0, 0, 0.5), 0.005)     # comprimento nulo
    with pytest.raises(ValueError):
        Condutor((0, 0, -0.1), (1, 0, 0.5), 0.005)    # z < 0 (acima da superficie)
    with pytest.raises(ValueError):
        Condutor((0, 0, 0.5), (1, 0, 0.5), 0.0)       # raio nao positivo


def test_segmentar_conta_e_versores():
    el = Eletrodo([Condutor((0, 0, 0.5), (10, 0, 0.5), 0.005)])
    s = el.segmentar(comp_alvo=1.0)
    assert s.n == 10
    assert np.allclose(2.0 * s.meia, 1.0)                       # cada seg 1 m
    assert np.allclose(np.linalg.norm(s.dir, axis=1), 1.0)     # versores
    assert (2.0 * s.meia).sum() == pytest.approx(10.0)         # comprimento total
    assert s.mid[0] == pytest.approx([0.5, 0.0, 0.5])          # 1o ponto medio


def test_segmentar_comp_alvo_nao_divisor():
    el = Eletrodo([Condutor((0, 0, 0.5), (10, 0, 0.5), 0.005)])
    s = el.segmentar(comp_alvo=3.0)
    assert s.n == 4                                            # ceil(10/3)
    assert np.allclose(2.0 * s.meia, 2.5)                      # 10/4


def test_segmentar_quebra_em_interface():
    # haste vertical 0 -> 4 m, interface em z = 2: nenhum segmento cruza z=2.
    el = Eletrodo([Condutor((0, 0, 0.0), (0, 0, 4.0), 0.005)])
    s = el.segmentar(comp_alvo=1.0, interfaces=[2.0])
    a = s.mid - s.meia[:, None] * s.dir
    b = s.mid + s.meia[:, None] * s.dir
    zlo = np.minimum(a[:, 2], b[:, 2])
    zhi = np.maximum(a[:, 2], b[:, 2])
    assert not np.any((zlo < 2.0 - 1e-9) & (zhi > 2.0 + 1e-9))
    # camada coerente com a profundidade do ponto medio
    assert np.all(s.camada[s.mid[:, 2] < 2.0] == 0)
    assert np.all(s.camada[s.mid[:, 2] > 2.0] == 1)


def test_malha_retangular_comprimento_condutores():
    el = Eletrodo.malha_retangular(comprimento_x=70.0, comprimento_y=70.0,
                                   espac_x=7.0, espac_y=7.0, prof_h=0.5, d=0.01)
    # 11 linhas em x + 11 em y, cada 70 m -> Lc = 1540 m.
    Lc = sum(c.comprimento for c in el.condutores)
    assert Lc == pytest.approx(1540.0)
    assert all(c.p1[2] == pytest.approx(0.5) for c in el.condutores)  # na prof_h


def test_malha_retangular_hastes():
    el = Eletrodo.malha_retangular(70.0, 70.0, 7.0, 7.0, 0.5, 0.01,
                                   n_hastes=4, comp_haste=7.5)
    total = sum(c.comprimento for c in el.condutores)
    assert total == pytest.approx(1540.0 + 4 * 7.5)           # grade + hastes
    # hastes sao verticais e descem a partir da malha (z cresce)
    verticais = [c for c in el.condutores
                 if c.p1[0] == c.p2[0] and c.p1[1] == c.p2[1]]
    assert len(verticais) == 4
    assert all(max(c.p1[2], c.p2[2]) == pytest.approx(0.5 + 7.5) for c in verticais)


def test_malha_retangular_hastes_demais_que_nos():
    # 10x10 com espac. 5 -> 8 nos de perimetro; 12 hastes coincidiriam.
    with pytest.raises(ValueError):
        Eletrodo.malha_retangular(10.0, 10.0, 5.0, 5.0, 0.5, 0.01,
                                  n_hastes=12, comp_haste=2.4)


def test_malha_retangular_hastes_em_posicoes_unicas():
    el = Eletrodo.malha_retangular(10.0, 10.0, 5.0, 5.0, 0.5, 0.01,
                                   n_hastes=8, comp_haste=2.4)
    hastes = [c for c in el.condutores
              if c.p1[0] == c.p2[0] and c.p1[1] == c.p2[1]]
    assert len(hastes) == 8
    assert len({(c.p1[0], c.p1[1]) for c in hastes}) == 8
