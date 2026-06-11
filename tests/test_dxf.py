"""Testes do importador de DXF (earthsolver.dxf)."""

import ezdxf
import pytest

from earthsolver.dxf import escanear_layers, from_dxf, wizard_mapa
from earthsolver.numerico import Eletrodo


def _criar_dxf(caminho, insunits=6):
    """Cria um DXF de teste: malha de 2 linhas (layer MALHA) + 2 hastes
    (POINT na layer HASTES). insunits no cabecalho (6 = metros)."""
    doc = ezdxf.new()
    doc.header["$INSUNITS"] = insunits
    doc.layers.add("MALHA")
    doc.layers.add("HASTES")
    msp = doc.modelspace()
    msp.add_line((0, 0), (10, 0), dxfattribs={"layer": "MALHA"})
    msp.add_line((0, 10), (10, 10), dxfattribs={"layer": "MALHA"})
    msp.add_point((0, 0), dxfattribs={"layer": "HASTES"})
    msp.add_point((10, 0), dxfattribs={"layer": "HASTES"})
    doc.saveas(caminho)
    return caminho


MAPA = {
    "padrao": {"prof": 0.5, "raio": 0.005},
    "layers": {
        "MALHA": {"prof": 0.5, "raio": 0.005},
        "HASTES": {"rod": True, "prof": 0.5, "comp": 3.0, "raio": 0.008},
    },
}


def test_linhas_viram_condutores_na_profundidade(tmp_path):
    arq = _criar_dxf(tmp_path / "m.dxf")
    el = from_dxf(arq, mapa=MAPA)
    assert isinstance(el, Eletrodo)
    linhas = [c for c in el.condutores if c.p1[2] == c.p2[2]]
    assert len(linhas) == 2
    c0 = linhas[0]
    assert c0.p1 == (0.0, 0.0, 0.5)
    assert c0.p2 == (10.0, 0.0, 0.5)
    assert c0.raio == 0.005


def test_hastes_viram_condutores_verticais(tmp_path):
    arq = _criar_dxf(tmp_path / "m.dxf")
    el = from_dxf(arq, mapa=MAPA)
    hastes = [c for c in el.condutores if c.p1[2] != c.p2[2]]
    assert len(hastes) == 2
    h = hastes[0]
    assert h.p1 == (0.0, 0.0, 0.5)
    assert h.p2 == (0.0, 0.0, 3.5)        # prof + comp
    assert h.raio == 0.008


def test_escala_converte_unidades(tmp_path):
    arq = _criar_dxf(tmp_path / "m.dxf")
    el = from_dxf(arq, mapa=MAPA, escala=0.001)   # desenho em mm
    linha = next(c for c in el.condutores if c.p1[2] == c.p2[2])
    assert linha.p2 == pytest.approx((0.01, 0.0, 0.5))   # 10 mm -> 0.01 m
    assert linha.p1[2] == 0.5                      # profundidade NAO escala


def test_layer_ausente_usa_padrao(tmp_path):
    arq = _criar_dxf(tmp_path / "m.dxf")
    mapa = {"padrao": {"prof": 1.2, "raio": 0.01}, "layers": {}}
    el = from_dxf(arq, mapa=mapa)
    linha = next(c for c in el.condutores if c.p1[2] == c.p2[2])
    assert linha.p1[2] == 1.2
    assert linha.raio == 0.01


def test_lwpolyline_gera_um_condutor_por_segmento(tmp_path):
    doc = ezdxf.new()
    doc.layers.add("MALHA")
    msp = doc.modelspace()
    msp.add_lwpolyline([(0, 0), (10, 0), (10, 10)], dxfattribs={"layer": "MALHA"})
    arq = tmp_path / "p.dxf"
    doc.saveas(arq)
    el = from_dxf(arq, mapa={"padrao": {"prof": 0.5, "raio": 0.005}, "layers": {}})
    assert len(el.condutores) == 2                 # 3 vertices -> 2 segmentos


def test_escanear_layers_lista_tipos(tmp_path):
    arq = _criar_dxf(tmp_path / "m.dxf")
    doc = ezdxf.readfile(arq)
    scan = escanear_layers(doc)
    assert scan["MALHA"]["tipos"] == {"LINE"}
    assert scan["MALHA"]["n"] == 2
    assert scan["HASTES"]["tipos"] == {"POINT"}
    assert scan["HASTES"]["n"] == 2


def _ler_roteirizado(respostas):
    it = iter(respostas)
    return lambda _prompt="": next(it)


def test_wizard_sugere_haste_para_layer_de_pontos(tmp_path):
    arq = _criar_dxf(tmp_path / "m.dxf")
    # escala (blank->sug), MALHA: tipo blank(->condutor), prof blank, raio blank,
    # HASTES: tipo blank(->haste), prof blank, raio blank, comp blank
    ler = _ler_roteirizado(["", "", "", "", "", "", "", ""])
    mapa = wizard_mapa(arq, ler=ler, escrever=lambda *a, **k: None)
    assert mapa["layers"]["HASTES"]["rod"] is True
    assert mapa["layers"]["MALHA"].get("rod", False) is False


def test_wizard_le_insunits_para_sugerir_escala(tmp_path):
    arq = _criar_dxf(tmp_path / "m.dxf", insunits=4)   # 4 = milimetros
    ler = _ler_roteirizado(["", "", "", "", "", "", "", ""])
    mapa = wizard_mapa(arq, ler=ler, escrever=lambda *a, **k: None)
    assert mapa["escala"] == pytest.approx(0.001)


def test_from_dxf_honra_escala_do_mapa(tmp_path):
    arq = _criar_dxf(tmp_path / "m.dxf")
    el = from_dxf(arq, mapa={**MAPA, "escala": 0.001})
    linha = next(c for c in el.condutores if c.p1[2] == c.p2[2])
    assert linha.p2 == pytest.approx((0.01, 0.0, 0.5))   # 10 mm -> 0.01 m


def test_escala_explicita_sobrepoe_a_do_mapa(tmp_path):
    arq = _criar_dxf(tmp_path / "m.dxf")
    el = from_dxf(arq, mapa={**MAPA, "escala": 0.001}, escala=1.0)
    linha = next(c for c in el.condutores if c.p1[2] == c.p2[2])
    assert linha.p2 == pytest.approx((10.0, 0.0, 0.5))


def test_wizard_repergunta_resposta_invalida(tmp_path):
    arq = _criar_dxf(tmp_path / "m.dxf")
    # escala, MALHA: "x" (invalida) -> repergunta -> "c", prof, raio;
    # HASTES: "h", prof, raio, comp
    ler = _ler_roteirizado(["1", "x", "c", "0.7", "0.006",
                            "h", "0.7", "0.01", "2.4"])
    mapa = wizard_mapa(arq, ler=ler, escrever=lambda *a, **k: None)
    assert mapa["layers"]["MALHA"] == {"prof": 0.7, "raio": 0.006}


def test_wizard_respostas_explicitas(tmp_path):
    arq = _criar_dxf(tmp_path / "m.dxf")
    # escala=1, MALHA condutor prof 0.7 raio 0.006, HASTES haste prof 0.7 raio 0.01 comp 2.4
    ler = _ler_roteirizado(["1", "c", "0.7", "0.006", "h", "0.7", "0.01", "2.4"])
    mapa = wizard_mapa(arq, ler=ler, escrever=lambda *a, **k: None)
    assert mapa["layers"]["MALHA"] == {"prof": 0.7, "raio": 0.006}
    assert mapa["layers"]["HASTES"] == {"rod": True, "prof": 0.7, "raio": 0.01, "comp": 2.4}
