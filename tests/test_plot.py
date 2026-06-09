"""Smoke tests dos graficos (earthsolver.plot). Usa o backend Agg via Figure."""

import os

import numpy as np

from earthsolver.numerico import Condutor, Eletrodo
from earthsolver.plot import plot_malha, plot_potencial, salvar


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
