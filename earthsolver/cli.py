"""Interface de linha de comando do earthsolver (argparse, stdlib).

Subcomandos:
  estratificar  - estratifica o solo a partir de um CSV de sondagem de Wenner.
  malha         - simula uma malha dada um solo e uma geometria (JSON).
  analisar      - pipeline completo a partir de um arquivo de projeto JSON.
"""

import argparse
import csv
import json
import sys

from .estratificacao import Estratificador
from .malha import EstudoAterramento, Malha
from .numerico import Condutor, Eletrodo, EstudoNumerico
from .solo import ModeloSolo


def _ler_sondagem_csv(caminho):
    """Le um CSV de sondagem. Colunas aceitas: spacing,resistance ou spacing,rho_a.

    Aceita cabecalho. A primeira coluna e o espacamento `a`; a segunda e a
    resistencia (default) ou a resistividade aparente se o cabecalho indicar.
    """
    espac, valores, col_rho = [], [], False
    with open(caminho, newline="") as f:
        leitor = csv.reader(f)
        linhas = [linha for linha in leitor if linha and linha[0].strip()]
    if not linhas:
        raise ValueError(f"arquivo de sondagem vazio: {caminho}")
    inicio = 0
    cab = [c.strip().lower() for c in linhas[0]]
    try:
        float(linhas[0][0])
    except ValueError:
        inicio = 1  # primeira linha e cabecalho
        col_rho = any("rho" in c or "resistivid" in c for c in cab)
    for linha in linhas[inicio:]:
        espac.append(float(linha[0]))
        valores.append(float(linha[1]))
    return espac, valores, col_rho


def _estratificar(args):
    espac, valores, col_rho = _ler_sondagem_csv(args.entrada)
    if col_rho:
        estrat = Estratificador(espac, resistividades=valores)
    else:
        estrat = Estratificador(espac, resistencias=valores)
    if args.camadas == "auto":
        estrat.auto_estratificar(max_camadas=args.max_camadas)
    else:
        estrat.estratificar(int(args.camadas))
    estrat.imprimir_modelo()
    if args.saida:
        estrat.exportar(args.saida)
    return estrat


def _carregar_solo(caminho):
    with open(caminho) as f:
        return ModeloSolo.from_dict(json.load(f))


def _carregar_malha(caminho):
    with open(caminho) as f:
        d = json.load(f)
    return Malha(**d)


def _malha(args):
    solo = _carregar_solo(args.solo)
    malha = _carregar_malha(args.malha)
    estudo = EstudoAterramento(solo, malha, Ig=args.ig, t=args.t, peso=args.peso)
    estudo.resolver()
    estudo.imprimir_resultado()
    if args.exportar:
        estudo.exportar(args.exportar)
    return estudo


def _carregar_eletrodo(caminho):
    """Le a geometria do eletrodo. Aceita {'condutores':[...]} ou {'malha_retangular':{...}}."""
    with open(caminho) as f:
        d = json.load(f)
    if "malha_retangular" in d:
        return Eletrodo.malha_retangular(**d["malha_retangular"])
    if "condutores" in d:
        return Eletrodo([Condutor(c["p1"], c["p2"], c["raio"])
                         for c in d["condutores"]])
    raise ValueError("eletrodo deve conter 'condutores' ou 'malha_retangular'")


def _numerico(args):
    solo = _carregar_solo(args.solo)
    eletrodo = _carregar_eletrodo(args.eletrodo)
    estudo = EstudoNumerico(solo, eletrodo, Ig=args.ig, t=args.t, peso=args.peso,
                            comp_alvo=args.comp_alvo, rho_s=args.rho_s, h_s=args.h_s)
    estudo.resolver()
    estudo.imprimir_resultado()
    if args.exportar:
        estudo.exportar(args.exportar)
    if args.raster:
        estudo.exportar_raster(args.raster)
    return estudo


def _analisar(args):
    """Projeto JSON: {solo:{...} ou sondagem:{...}, malha:{...}, falta:{...}}."""
    with open(args.projeto) as f:
        proj = json.load(f)

    if "solo" in proj:
        solo = ModeloSolo.from_dict(proj["solo"])
    elif "sondagem" in proj:
        s = proj["sondagem"]
        estrat = Estratificador(
            s["espacamentos"],
            resistencias=s.get("resistencias"),
            resistividades=s.get("resistividades"),
        )
        camadas = s.get("camadas", "auto")
        if camadas == "auto":
            estrat.auto_estratificar(max_camadas=s.get("max_camadas", 4))
        else:
            estrat.estratificar(int(camadas))
        estrat.imprimir_modelo()
        solo = estrat.modelo
    else:
        raise ValueError("projeto deve conter 'solo' ou 'sondagem'")

    malha = Malha(**proj["malha"])
    falta = proj.get("falta", {})
    estudo = EstudoAterramento(
        solo, malha,
        Ig=falta.get("Ig", proj.get("Ig")),
        t=falta.get("t", proj.get("t")),
        peso=falta.get("peso", proj.get("peso", 70)),
    )
    estudo.resolver()
    estudo.imprimir_resultado()
    if args.exportar:
        estudo.exportar(args.exportar)
    return estudo


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="earthsolver",
        description="Suite de analise de aterramento (IEEE Std 80 / NBR).",
    )
    sub = parser.add_subparsers(dest="comando", required=True)

    p_e = sub.add_parser("estratificar", help="estratifica o solo de um CSV de Wenner")
    p_e.add_argument("--entrada", required=True, help="CSV: spacing,resistance|rho_a")
    p_e.add_argument("--camadas", default="auto", help="'auto' ou numero de camadas")
    p_e.add_argument("--max-camadas", type=int, default=4, dest="max_camadas")
    p_e.add_argument("--saida", help="arquivo JSON de saida do modelo de solo")
    p_e.set_defaults(func=_estratificar)

    p_m = sub.add_parser("malha", help="simula uma malha dado solo + geometria")
    p_m.add_argument("--solo", required=True, help="JSON do modelo de solo")
    p_m.add_argument("--malha", required=True, help="JSON da geometria da malha")
    p_m.add_argument("--ig", type=float, required=True, help="corrente de malha (A)")
    p_m.add_argument("--t", type=float, required=True, help="duracao da falta (s)")
    p_m.add_argument("--peso", type=int, default=70, help="peso corporeo 50 ou 70")
    p_m.add_argument("--exportar", help="arquivo JSON de saida")
    p_m.set_defaults(func=_malha)

    p_a = sub.add_parser("analisar", help="pipeline completo a partir de um projeto JSON")
    p_a.add_argument("projeto", help="arquivo de projeto JSON")
    p_a.add_argument("--exportar", help="arquivo JSON de saida")
    p_a.set_defaults(func=_analisar)

    p_n = sub.add_parser("numerico",
                         help="solver numerico de segmentacao (geometria explicita)")
    p_n.add_argument("--eletrodo", required=True,
                     help="JSON: {condutores:[...]} ou {malha_retangular:{...}}")
    p_n.add_argument("--solo", required=True, help="JSON do modelo de solo")
    p_n.add_argument("--ig", type=float, required=True, help="corrente de malha (A)")
    p_n.add_argument("--t", type=float, required=True, help="duracao da falta (s)")
    p_n.add_argument("--peso", type=int, default=70, help="peso corporeo 50 ou 70")
    p_n.add_argument("--comp-alvo", type=float, default=2.0, dest="comp_alvo",
                     help="comprimento alvo de segmento (m)")
    p_n.add_argument("--rho-s", type=float, default=None, dest="rho_s",
                     help="resistividade da camada superficial (brita), Ohm.m")
    p_n.add_argument("--h-s", type=float, default=0.1, dest="h_s",
                     help="espessura da camada superficial (m)")
    p_n.add_argument("--exportar", help="arquivo JSON de saida do resultado")
    p_n.add_argument("--raster", help="arquivo JSON do mapa de potencial de superficie")
    p_n.set_defaults(func=_numerico)

    args = parser.parse_args(argv)
    try:
        args.func(args)
    except (ValueError, OSError, KeyError) as exc:
        print(f"erro: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
