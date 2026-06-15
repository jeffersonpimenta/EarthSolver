"""Testes do subcomando CLI 'numerico'."""

import json

from earthsolver.cli import main


def _escrever(p, obj):
    p.write_text(json.dumps(obj))


_MALHA = {"malha_retangular": {"comprimento_x": 70, "comprimento_y": 70,
                               "espac_x": 7, "espac_y": 7, "prof_h": 0.5, "d": 0.01}}


def test_cli_numerico_malha_retangular(tmp_path):
    solo = tmp_path / "solo.json"
    _escrever(solo, {"rho": [400.0], "espessura": []})
    cond = tmp_path / "cond.json"
    _escrever(cond, _MALHA)
    out = tmp_path / "out.json"
    rc = main(["numerico", "--eletrodo", str(cond), "--solo", str(solo),
               "--ig", "1908", "--t", "0.5", "--comp-alvo", "7",
               "--exportar", str(out)])
    assert rc == 0
    dados = json.loads(out.read_text())
    assert dados["resultado"]["Rg"] > 0
    assert "aprovado" in dados["resultado"]


def test_cli_numerico_condutores_explicitos(tmp_path):
    solo = tmp_path / "solo.json"
    _escrever(solo, {"rho": [100.0], "espessura": []})
    cond = tmp_path / "cond.json"
    _escrever(cond, {"condutores": [{"p1": [0, 0, 0.0], "p2": [0, 0, 3.0],
                                     "raio": 0.01}]})
    rc = main(["numerico", "--eletrodo", str(cond), "--solo", str(solo),
               "--ig", "10", "--t", "0.5", "--comp-alvo", "0.5"])
    assert rc == 0


def test_cli_numerico_exporta_raster(tmp_path):
    solo = tmp_path / "solo.json"
    _escrever(solo, {"rho": [400.0], "espessura": []})
    cond = tmp_path / "cond.json"
    _escrever(cond, _MALHA)
    rast = tmp_path / "pot.json"
    rc = main(["numerico", "--eletrodo", str(cond), "--solo", str(solo),
               "--ig", "1908", "--t", "0.5", "--comp-alvo", "7",
               "--raster", str(rast)])
    assert rc == 0 and rast.exists()
    assert "phi" in json.loads(rast.read_text())


def test_cli_numerico_plots_seguranca(tmp_path):
    solo = tmp_path / "solo.json"
    _escrever(solo, {"rho": [400.0], "espessura": []})
    cond = tmp_path / "cond.json"
    _escrever(cond, _MALHA)
    toque = tmp_path / "toque.png"
    passo = tmp_path / "passo.png"
    margem = tmp_path / "margem.png"
    perfis = tmp_path / "perfis.png"
    corrente = tmp_path / "corrente.png"
    rc = main(["numerico", "--eletrodo", str(cond), "--solo", str(solo),
               "--ig", "1908", "--t", "0.5", "--comp-alvo", "7",
               "--plot-toque", str(toque), "--plot-passo", str(passo),
               "--plot-margem", str(margem), "--plot-perfis", str(perfis),
               "--plot-corrente", str(corrente)])
    assert rc == 0
    for png in (toque, passo, margem, perfis, corrente):
        assert png.exists() and png.stat().st_size > 0


def test_cli_numerico_plots_3d_e_convergencia(tmp_path):
    solo = tmp_path / "solo.json"
    _escrever(solo, {"rho": [400.0], "espessura": []})
    cond = tmp_path / "cond.json"
    _escrever(cond, _MALHA)
    toque3d = tmp_path / "toque3d.png"
    passo3d = tmp_path / "passo3d.png"
    conv = tmp_path / "conv.png"
    rc = main(["numerico", "--eletrodo", str(cond), "--solo", str(solo),
               "--ig", "1908", "--t", "0.5", "--comp-alvo", "7",
               "--plot-toque-3d", str(toque3d), "--plot-passo-3d", str(passo3d),
               "--plot-convergencia", str(conv), "--conv-alvos", "14,7"])
    assert rc == 0
    for png in (toque3d, passo3d, conv):
        assert png.exists() and png.stat().st_size > 0
