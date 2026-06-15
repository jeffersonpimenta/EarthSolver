"""Testes do campo de potencial de superficie, toque/passo e raster (numerico)."""

import numpy as np

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


def test_raster_toque_passo_existem_e_consistentes():
    est = _est()
    est.resolver()
    X, Y, _ = est.raster
    Xt, Yt, toque = est.raster_toque
    Xs, Ys, passo = est.raster_passo
    # mesma grade do raster de potencial
    assert toque.shape == passo.shape == X.shape
    assert np.array_equal(Xt, X) and np.array_equal(Yt, Y)
    assert np.array_equal(Xs, X) and np.array_equal(Ys, Y)
    assert np.all(np.isfinite(toque)) and np.all(np.isfinite(passo))


def test_campos_reproduzem_escalares_em_es():
    est = _est()
    est.resolver()
    Xt, Yt, toque = est.raster_toque
    _, _, passo = est.raster_passo
    # Es = pior caso do campo de passo (toda a area)
    assert passo.max() == est.Es
    # Em = pior toque dentro da projecao da malha (bbox dos segmentos)
    pts = np.vstack([est._A[:, :2], est._B[:, :2]])
    xmin, ymin = pts.min(axis=0)
    xmax, ymax = pts.max(axis=0)
    dentro = (Xt >= xmin) & (Xt <= xmax) & (Yt >= ymin) & (Yt <= ymax)
    assert np.isclose(toque[dentro].max(), est.Em)


def test_dados_corrente_devolve_A_B_I():
    est = _est()
    est.resolver()
    A, B, I = est.dados_corrente()
    M = est.segs.n
    assert A.shape == (M, 3) and B.shape == (M, 3)
    assert I.shape == (M,)
    # corrente total drenada = corrente injetada Ig
    assert np.isclose(I.sum(), est.Ig)
