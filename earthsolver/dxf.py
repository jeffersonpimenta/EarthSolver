"""Importacao de geometria de aterramento a partir de arquivos DXF.

Um DXF de projeto e plano (so x, y). A profundidade de enterramento, o raio do
condutor e quais layers sao hastes verticais vem de um *mapa de layers* (dict):

    {
      "padrao": {"prof": 0.5, "raio": 0.005},
      "layers": {
        "MALHA":  {"prof": 0.5, "raio": 0.005},
        "HASTES": {"rod": true, "prof": 0.5, "comp": 7.5, "raio": 0.008}
      }
    }

`from_dxf` converte para um `Eletrodo` (reaproveita `Condutor`/`Eletrodo`).
Quando o usuario nao tem um mapa, `wizard_mapa` escaneia as layers e o monta
interativamente. `ezdxf` e importado de forma preguicosa (so aqui).
"""

import math
import warnings

from .numerico import Condutor, Eletrodo

LINHA_TIPOS = {"LINE", "LWPOLYLINE", "POLYLINE"}
PONTO_TIPOS = {"POINT", "INSERT", "CIRCLE"}

# Codigo $INSUNITS do DXF -> fator para metros (1 desenho = ? m).
_INSUNITS_ESCALA = {1: 0.0254, 2: 0.3048, 4: 0.001, 5: 0.01, 6: 1.0}

_PADRAO = {"prof": 0.5, "raio": 0.005}


def escanear_layers(doc) -> dict:
    """Varre o modelspace e devolve {layer: {'tipos': set, 'n': contagem}}.

    A ordem das chaves segue a ordem de encontro das entidades no desenho.
    """
    scan = {}
    for e in doc.modelspace():
        d = scan.setdefault(e.dxf.layer, {"tipos": set(), "n": 0})
        d["tipos"].add(e.dxftype())
        d["n"] += 1
    return scan


def from_dxf(caminho, mapa=None, padrao=None, escala=1.0) -> Eletrodo:
    """Le um DXF e devolve o `Eletrodo` correspondente.

    mapa: dict {'padrao':..., 'layers':{...}} (ou None -> tudo com `padrao`).
    padrao: config da layer nao listada, se `mapa` nao trouxer 'padrao'.
    escala: multiplica as coordenadas x, y (CAD em mm -> 0.001). prof/raio/comp
    sao tomados em metros (nao escalam).
    """
    import ezdxf

    doc = ezdxf.readfile(str(caminho))
    layers = (mapa or {}).get("layers", {})
    padrao_cfg = (mapa or {}).get("padrao") or padrao or _PADRAO

    condutores = []
    for e in doc.modelspace():
        tipo = e.dxftype()
        cfg = layers.get(e.dxf.layer, padrao_cfg)
        prof = float(cfg.get("prof", _PADRAO["prof"]))
        raio = float(cfg.get("raio", _PADRAO["raio"]))

        if cfg.get("rod"):
            if tipo in PONTO_TIPOS:
                comp = float(cfg.get("comp", 0.0))
                if comp <= 0:
                    raise ValueError(
                        f"layer de haste '{e.dxf.layer}' precisa de 'comp' > 0")
                x, y = _xy_ponto(e, tipo)
                condutores.append(Condutor((x * escala, y * escala, prof),
                                           (x * escala, y * escala, prof + comp),
                                           raio))
            elif tipo in LINHA_TIPOS:
                warnings.warn(
                    f"{tipo} em layer de haste '{e.dxf.layer}' ignorada")
        else:
            if tipo == "LINE":
                (x1, y1), (x2, y2) = _xy_linha(e)
                _add_seg(condutores, x1, y1, x2, y2, prof, raio, escala)
            elif tipo in ("LWPOLYLINE", "POLYLINE"):
                pts, fechado = _pontos_polilinha(e)
                segs = list(zip(pts, pts[1:]))
                if fechado and len(pts) > 2:
                    segs.append((pts[-1], pts[0]))
                for (x1, y1), (x2, y2) in segs:
                    _add_seg(condutores, x1, y1, x2, y2, prof, raio, escala)
            elif tipo in PONTO_TIPOS:
                warnings.warn(
                    f"{tipo} em '{e.dxf.layer}' ignorado "
                    f"(marque a layer como haste para virar haste)")
            elif tipo in ("ARC", "SPLINE", "ELLIPSE"):
                warnings.warn(
                    f"{tipo} em '{e.dxf.layer}' ignorada (fora do escopo v1)")

    return Eletrodo(condutores)


def wizard_mapa(caminho, ler=input, escrever=print) -> dict:
    """Escaneia as layers do DXF e monta o mapa interativamente.

    ler/escrever sao injetaveis (default input/print) para teste sem TTY.
    Devolve {'escala':.., 'padrao':{...}, 'layers':{...}}.
    """
    import ezdxf

    doc = ezdxf.readfile(str(caminho))
    scan = escanear_layers(doc)

    escrever(f"Assistente de mapa de layers - {caminho}")
    insunits = int(doc.header.get("$INSUNITS", 0) or 0)
    sug_esc = _INSUNITS_ESCALA.get(insunits, 1.0)
    escala = _ler_float(ler, f"Escala p/ metros (1 desenho = ? m) [{sug_esc}]: ",
                        sug_esc)

    layers = {}
    for lay, info in scan.items():
        tipos = info["tipos"]
        sug = "h" if tipos <= PONTO_TIPOS else "c"
        escrever(f"Layer '{lay}': {sorted(tipos)} ({info['n']} entidades)")
        resp = (ler(f"  [c]ondutor / [h]aste / [i]gnorar [{sug}]: ") or sug)
        resp = resp.strip().lower()[:1] or sug
        if resp == "i":
            continue
        prof = _ler_float(ler, "  profundidade (m) [0.5]: ", 0.5)
        raio = _ler_float(ler, "  raio do condutor (m) [0.005]: ", 0.005)
        if resp == "h":
            comp = _ler_float(ler, "  comprimento da haste (m) [2.4]: ", 2.4)
            layers[lay] = {"rod": True, "prof": prof, "raio": raio, "comp": comp}
        else:
            layers[lay] = {"prof": prof, "raio": raio}

    return {"escala": escala, "padrao": dict(_PADRAO), "layers": layers}


# --------------------------------------------------------------- helpers
def _ler_float(ler, prompt, default):
    resp = (ler(prompt) or "").strip()
    return float(resp) if resp else float(default)


def _add_seg(lst, x1, y1, x2, y2, prof, raio, escala):
    p1 = (x1 * escala, y1 * escala, prof)
    p2 = (x2 * escala, y2 * escala, prof)
    if math.dist(p1, p2) > 1e-12:
        lst.append(Condutor(p1, p2, raio))


def _xy_linha(e):
    s, t = e.dxf.start, e.dxf.end
    return (float(s.x), float(s.y)), (float(t.x), float(t.y))


def _xy_ponto(e, tipo):
    v = {"POINT": "location", "INSERT": "insert", "CIRCLE": "center"}[tipo]
    p = getattr(e.dxf, v)
    return float(p.x), float(p.y)


def _pontos_polilinha(e):
    if e.dxftype() == "LWPOLYLINE":
        pts = [(float(x), float(y)) for x, y in e.get_points("xy")]
        return pts, bool(e.closed)
    pts = [(float(v.dxf.location.x), float(v.dxf.location.y)) for v in e.vertices]
    return pts, bool(e.is_closed)
