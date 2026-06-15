"""Smoke tests dos graficos (earthsolver.plot). Usa o backend Agg via Figure."""


import numpy as np

from earthsolver.numerico import Condutor, Eletrodo
from earthsolver.plot import (
    plot_corrente,
    plot_malha,
    plot_malha_3d,
    plot_margem,
    plot_potencial,
    plot_potencial_3d,
    plot_perfis,
    plot_tensao_passo,
    plot_tensao_toque,
    salvar,
)


def _eletrodo():
    return Eletrodo([
        Condutor((0, 0, 0.5), (10, 0, 0.5), 0.005),    # condutor horizontal
        Condutor((0, 0, 0.5), (0, 0, 3.0), 0.008),     # haste vertical
    ])


def _raster():
    xs = np.linspace(0, 10, 11)
    ys = np.linspace(0, 10, 11)
    X, Y = np.meshgrid(xs, ys)
    Phi = 1000.0 * np.exp(-((X - 5) ** 2 + (Y - 5) ** 2) / 10.0)
    return X, Y, Phi


def _raster_toque():
    X, Y, Phi = _raster()
    return X, Y, 1000.0 - Phi            # GPR=1000 -> toque = GPR - Phi


def _raster_passo():
    X, Y, _ = _raster()
    passo = 300.0 * np.exp(-((X - 5) ** 2 + (Y - 5) ** 2) / 8.0)
    return X, Y, passo


def _segmentos():
    A = np.array([[0.0, 0.0, 0.5], [5.0, 0.0, 0.5], [0.0, 0.0, 0.5]])
    B = np.array([[5.0, 0.0, 0.5], [10.0, 0.0, 0.5], [0.0, 0.0, 3.0]])
    I = np.array([0.8, 0.5, 0.2])
    return A, B, I


def test_plot_malha_salva_png(tmp_path):
    png = tmp_path / "malha.png"
    ax = plot_malha(_eletrodo())
    salvar(ax, png)
    assert png.exists() and png.stat().st_size > 0


def test_plot_potencial_salva_png(tmp_path):
    png = tmp_path / "pot.png"
    ax = plot_potencial(_raster(), eletrodo=_eletrodo())
    salvar(ax, png)
    assert png.exists() and png.stat().st_size > 0


def test_plot_malha_3d_salva_png(tmp_path):
    png = tmp_path / "malha_3d.png"
    ax = plot_malha_3d(_eletrodo())
    salvar(ax, png)
    assert png.exists() and png.stat().st_size > 0


def test_plot_potencial_3d_salva_png(tmp_path):
    png = tmp_path / "pot_3d.png"
    ax = plot_potencial_3d(_raster(), eletrodo=_eletrodo())
    salvar(ax, png)
    assert png.exists() and png.stat().st_size > 0


def test_plot_tensao_toque_salva_png(tmp_path):
    png = tmp_path / "toque.png"
    ax = plot_tensao_toque(_raster_toque(), limite=500.0, eletrodo=_eletrodo())
    salvar(ax, png)
    assert png.exists() and png.stat().st_size > 0


def test_plot_tensao_toque_sem_limite(tmp_path):
    png = tmp_path / "toque2.png"
    ax = plot_tensao_toque(_raster_toque())
    salvar(ax, png)
    assert png.exists() and png.stat().st_size > 0


def test_plot_tensao_passo_salva_png(tmp_path):
    png = tmp_path / "passo.png"
    ax = plot_tensao_passo(_raster_passo(), limite=200.0, eletrodo=_eletrodo())
    salvar(ax, png)
    assert png.exists() and png.stat().st_size > 0


def test_plot_margem_salva_png(tmp_path):
    png = tmp_path / "margem.png"
    ax = plot_margem(_raster_toque(), _raster_passo(),
                     E_toque=500.0, E_passo=200.0, eletrodo=_eletrodo())
    salvar(ax, png)
    assert png.exists() and png.stat().st_size > 0


def test_plot_perfis_salva_png(tmp_path):
    png = tmp_path / "perfis.png"
    obj = plot_perfis(_raster(), _raster_toque(), _raster_passo(),
                      E_toque=500.0, E_passo=200.0, GPR=1000.0)
    salvar(obj, png)
    assert png.exists() and png.stat().st_size > 0


def test_plot_corrente_salva_png(tmp_path):
    png = tmp_path / "corrente.png"
    A, B, I = _segmentos()
    ax = plot_corrente(A, B, I)
    salvar(ax, png)
    assert png.exists() and png.stat().st_size > 0
