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
    for num, linha in enumerate(linhas[inicio:], start=inicio + 1):
        if len(linha) < 2:
            raise ValueError(
                f"linha {num} do CSV deve ter 2 colunas "
                f"(espacamento, valor): {linha!r}")
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
                            comp_alvo=args.comp_alvo, rho_s=args.rho_s, h_s=args.h_s,
                            precisao=args.precisao, mem_mb=args.mem_mb)
    prog = None if getattr(args, "sem_progresso", False) else "barra"
    estudo.resolver(progresso=prog)
    estudo.imprimir_resultado()
    if args.exportar:
        estudo.exportar(args.exportar)
    if args.raster:
        estudo.exportar_raster(args.raster)
    if getattr(args, "plot", None):
        from . import plot
        plot.salvar(plot.plot_potencial(estudo.raster, eletrodo=estudo.eletrodo),
                    args.plot)
        print(f"Mapa de potencial exportado para {args.plot}.")
    if getattr(args, "plot_malha", None):
        from . import plot
        plot.salvar(plot.plot_malha(estudo.eletrodo), args.plot_malha)
        print(f"Vista da malha exportada para {args.plot_malha}.")
    if getattr(args, "plot_3d", None):
        from . import plot
        plot.salvar(plot.plot_potencial_3d(estudo.raster, eletrodo=estudo.eletrodo),
                    args.plot_3d)
        print(f"Elevacao do potencial exportada para {args.plot_3d}.")
    if getattr(args, "plot_malha_3d", None):
        from . import plot
        plot.salvar(plot.plot_malha_3d(estudo.eletrodo), args.plot_malha_3d)
        print(f"Perspectiva da malha exportada para {args.plot_malha_3d}.")
    r = estudo.resultado
    if getattr(args, "plot_toque", None):
        from . import plot
        plot.salvar(plot.plot_tensao_toque(estudo.raster_toque,
                    limite=r["E_toque"], eletrodo=estudo.eletrodo),
                    args.plot_toque)
        print(f"Mapa de tensao de toque exportado para {args.plot_toque}.")
    if getattr(args, "plot_passo", None):
        from . import plot
        plot.salvar(plot.plot_tensao_passo(estudo.raster_passo,
                    limite=r["E_passo"], eletrodo=estudo.eletrodo),
                    args.plot_passo)
        print(f"Mapa de tensao de passo exportado para {args.plot_passo}.")
    if getattr(args, "plot_margem", None):
        from . import plot
        plot.salvar(plot.plot_margem(estudo.raster_toque, estudo.raster_passo,
                    r["E_toque"], r["E_passo"], eletrodo=estudo.eletrodo),
                    args.plot_margem)
        print(f"Mapa de margem exportado para {args.plot_margem}.")
    if getattr(args, "plot_perfis", None):
        from . import plot
        plot.salvar(plot.plot_perfis(estudo.raster, estudo.raster_toque,
                    estudo.raster_passo, r["E_toque"], r["E_passo"], r["GPR"]),
                    args.plot_perfis)
        print(f"Perfis em corte exportados para {args.plot_perfis}.")
    if getattr(args, "plot_corrente", None):
        from . import plot
        A, B, I = estudo.dados_corrente()
        plot.salvar(plot.plot_corrente(A, B, I), args.plot_corrente)
        print(f"Distribuicao de corrente exportada para {args.plot_corrente}.")
    if getattr(args, "plot_toque_3d", None):
        from . import plot
        plot.salvar(plot.plot_campo_3d(estudo.raster_toque,
                    rotulo="Tensao de toque (V)", titulo="Tensao de toque (3D)",
                    eletrodo=estudo.eletrodo), args.plot_toque_3d)
        print(f"Superficie 3D de toque exportada para {args.plot_toque_3d}.")
    if getattr(args, "plot_passo_3d", None):
        from . import plot
        plot.salvar(plot.plot_campo_3d(estudo.raster_passo,
                    rotulo="Tensao de passo (V)", titulo="Tensao de passo (3D)",
                    eletrodo=estudo.eletrodo), args.plot_passo_3d)
        print(f"Superficie 3D de passo exportada para {args.plot_passo_3d}.")
    if getattr(args, "plot_convergencia", None):
        from . import plot
        from .numerico import estudo_convergencia
        if getattr(args, "conv_alvos", None):
            alvos = [float(v) for v in args.conv_alvos.split(",")]
        else:
            alvos = [args.comp_alvo * k for k in (3.0, 2.0, 1.0, 0.6)]
        dados = estudo_convergencia(solo, eletrodo, args.ig, args.t, alvos,
                                    progresso=prog,
                                    peso=args.peso, rho_s=args.rho_s, h_s=args.h_s,
                                    precisao=args.precisao, mem_mb=args.mem_mb)
        plot.salvar(plot.plot_convergencia(dados), args.plot_convergencia)
        print(f"Curva de convergencia exportada para {args.plot_convergencia}.")
    return estudo


def _dxf(args):
    """Converte um DXF em geometria de eletrodo ({condutores:[...]} JSON).

    Sem --mapa: roda o wizard (se a sessao for interativa e nao vier --sem-wizard);
    caso contrario aplica os defaults --prof/--raio a todas as layers.
    """
    from . import dxf, plot

    mapa = None
    if args.mapa:
        with open(args.mapa) as f:
            mapa = json.load(f)
    elif sys.stdin.isatty() and not args.sem_wizard:
        mapa = dxf.wizard_mapa(args.arquivo)
        destino = args.salvar_mapa
        if not destino:
            resp = input("Salvar este mapa para reusar? caminho [Enter pula]: ")
            destino = resp.strip() or None
        if destino:
            with open(destino, "w") as f:
                json.dump(mapa, f, indent=2)
            print(f"Mapa salvo em {destino}.")
    else:
        mapa = {"padrao": {"prof": args.prof, "raio": args.raio}, "layers": {}}

    eletrodo = dxf.from_dxf(args.arquivo, mapa=mapa, escala=args.escala)
    condutores = [{"p1": list(c.p1), "p2": list(c.p2), "raio": c.raio}
                  for c in eletrodo.condutores]
    with open(args.saida, "w") as f:
        json.dump({"condutores": condutores}, f, indent=2)
    print(f"{len(condutores)} condutores exportados para {args.saida}.")
    if args.plot_malha:
        plot.salvar(plot.plot_malha(eletrodo), args.plot_malha)
        print(f"Vista da malha exportada para {args.plot_malha}.")
    if args.plot_malha_3d:
        plot.salvar(plot.plot_malha_3d(eletrodo), args.plot_malha_3d)
        print(f"Perspectiva da malha exportada para {args.plot_malha_3d}.")
    return eletrodo


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
    Ig = falta.get("Ig", proj.get("Ig"))
    t = falta.get("t", proj.get("t"))
    if Ig is None or t is None:
        raise ValueError(
            "projeto deve definir 'Ig' e 't' (em 'falta' ou no topo)")
    estudo = EstudoAterramento(
        solo, malha, Ig=Ig, t=t,
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
    p_n.add_argument("--precisao", choices=("double", "single"), default="double",
                     help="precisao da montagem/solucao: double (float64) ou "
                          "single (float32, ~metade da RAM)")
    p_n.add_argument("--mem-mb", type=float, default=256.0, dest="mem_mb",
                     help="orcamento de RAM (MB) por bloco da montagem da matriz; "
                          "menor = menos memoria de pico (comp-alvo pequeno)")
    p_n.add_argument("--exportar", help="arquivo JSON de saida do resultado")
    p_n.add_argument("--raster", help="arquivo JSON do mapa de potencial de superficie")
    p_n.add_argument("--plot", help="PNG do mapa de potencial de superficie")
    p_n.add_argument("--plot-malha", dest="plot_malha", help="PNG da vista da malha")
    p_n.add_argument("--plot-3d", dest="plot_3d",
                     help="PNG da elevacao do potencial em 3D")
    p_n.add_argument("--plot-malha-3d", dest="plot_malha_3d",
                     help="PNG da perspectiva (3D) da malha")
    p_n.add_argument("--plot-toque", dest="plot_toque",
                     help="PNG do mapa de tensao de toque")
    p_n.add_argument("--plot-passo", dest="plot_passo",
                     help="PNG do mapa de tensao de passo")
    p_n.add_argument("--plot-margem", dest="plot_margem",
                     help="PNG do mapa de margem de seguranca (utilizacao %%)")
    p_n.add_argument("--plot-perfis", dest="plot_perfis",
                     help="PNG dos perfis em corte (potencial/toque/passo)")
    p_n.add_argument("--plot-corrente", dest="plot_corrente",
                     help="PNG da distribuicao de corrente por segmento")
    p_n.add_argument("--plot-toque-3d", dest="plot_toque_3d",
                     help="PNG da superficie 3D da tensao de toque")
    p_n.add_argument("--plot-passo-3d", dest="plot_passo_3d",
                     help="PNG da superficie 3D da tensao de passo")
    p_n.add_argument("--plot-convergencia", dest="plot_convergencia",
                     help="PNG da curva de convergencia (Rg/GPR vs nº segmentos)")
    p_n.add_argument("--conv-alvos", dest="conv_alvos",
                     help="comprimentos-alvo p/ a convergencia (ex.: 14,7,3.5,2)")
    p_n.add_argument("--sem-progresso", action="store_true", dest="sem_progresso",
                     help="desliga a barra de progresso/ETA no terminal")
    p_n.set_defaults(func=_numerico)

    p_d = sub.add_parser("dxf", help="converte um DXF em geometria de eletrodo")
    p_d.add_argument("arquivo", help="arquivo DXF de entrada")
    p_d.add_argument("--mapa", help="JSON do mapa de layers (sem ele, roda o wizard)")
    p_d.add_argument("--prof", type=float, default=0.5,
                     help="profundidade padrao (m) quando sem mapa")
    p_d.add_argument("--raio", type=float, default=0.005,
                     help="raio padrao (m) quando sem mapa")
    p_d.add_argument("--escala", type=float, default=None,
                     help="fator de escala p/ metros (desenho em mm -> 0.001)")
    p_d.add_argument("-o", "--saida", default="cond.json",
                     help="JSON de saida da geometria (default cond.json)")
    p_d.add_argument("--salvar-mapa", dest="salvar_mapa",
                     help="salva o mapa montado pelo wizard para reuso")
    p_d.add_argument("--sem-wizard", dest="sem_wizard", action="store_true",
                     help="nao roda o wizard; aplica --prof/--raio a tudo")
    p_d.add_argument("--plot-malha", dest="plot_malha",
                     help="PNG da vista da malha importada")
    p_d.add_argument("--plot-malha-3d", dest="plot_malha_3d",
                     help="PNG da perspectiva (3D) da malha importada")
    p_d.set_defaults(func=_dxf)

    args = parser.parse_args(argv)
    try:
        args.func(args)
    except (ValueError, OSError, KeyError) as exc:
        print(f"erro: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
