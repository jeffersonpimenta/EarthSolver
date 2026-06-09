"""Testes do subcomando CLI 'dxf' e das flags --plot do 'numerico'."""

import json

import ezdxf

from earthsolver.cli import main


def _criar_dxf(caminho):
    doc = ezdxf.new()
    doc.layers.add("MALHA")
    doc.layers.add("HASTES")
    msp = doc.modelspace()
    for y in (0, 7, 14):
        msp.add_line((0, y), (14, y), dxfattribs={"layer": "MALHA"})
    for x in (0, 7, 14):
        msp.add_line((x, 0), (x, 14), dxfattribs={"layer": "MALHA"})
    msp.add_point((0, 0), dxfattribs={"layer": "HASTES"})
    msp.add_point((14, 14), dxfattribs={"layer": "HASTES"})
    doc.saveas(caminho)
    return caminho


_MAPA = {
    "padrao": {"prof": 0.5, "raio": 0.005},
    "layers": {
        "MALHA": {"prof": 0.5, "raio": 0.005},
        "HASTES": {"rod": True, "prof": 0.5, "comp": 3.0, "raio": 0.008},
    },
}


def test_cli_dxf_com_mapa_gera_condutores_compativeis(tmp_path):
    dxf = _criar_dxf(tmp_path / "m.dxf")
    mapa = tmp_path / "layers.json"
    mapa.write_text(json.dumps(_MAPA))
    cond = tmp_path / "cond.json"
    rc = main(["dxf", str(dxf), "--mapa", str(mapa), "-o", str(cond)])
    assert rc == 0 and cond.exists()
    dados = json.loads(cond.read_text())
    assert len(dados["condutores"]) == 8        # 6 linhas + 2 hastes
    c = dados["condutores"][0]
    assert {"p1", "p2", "raio"} <= c.keys()

    # o cond.json gerado e diretamente consumivel pelo comando 'numerico'
    solo = tmp_path / "solo.json"
    solo.write_text(json.dumps({"rho": [400.0], "espessura": []}))
    rc2 = main(["numerico", "--eletrodo", str(cond), "--solo", str(solo),
                "--ig", "100", "--t", "0.5", "--comp-alvo", "7"])
    assert rc2 == 0


def test_cli_dxf_sem_wizard_usa_padrao(tmp_path):
    dxf = _criar_dxf(tmp_path / "m.dxf")
    cond = tmp_path / "cond.json"
    rc = main(["dxf", str(dxf), "--sem-wizard", "--prof", "0.8",
               "--raio", "0.01", "-o", str(cond)])
    assert rc == 0
    dados = json.loads(cond.read_text())
    assert all(c["p1"][2] == 0.8 for c in dados["condutores"])   # padrao aplicado


def test_cli_numerico_plot_gera_pngs(tmp_path):
    solo = tmp_path / "solo.json"
    solo.write_text(json.dumps({"rho": [400.0], "espessura": []}))
    cond = tmp_path / "cond.json"
    cond.write_text(json.dumps({"malha_retangular": {
        "comprimento_x": 70, "comprimento_y": 70, "espac_x": 7, "espac_y": 7,
        "prof_h": 0.5, "d": 0.01}}))
    pot = tmp_path / "pot.png"
    mal = tmp_path / "mal.png"
    rc = main(["numerico", "--eletrodo", str(cond), "--solo", str(solo),
               "--ig", "1908", "--t", "0.5", "--comp-alvo", "7",
               "--plot", str(pot), "--plot-malha", str(mal)])
    assert rc == 0
    assert pot.stat().st_size > 0 and mal.stat().st_size > 0
