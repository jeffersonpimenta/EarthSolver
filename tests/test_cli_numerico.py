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
