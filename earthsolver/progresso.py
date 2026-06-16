"""Relatorio de progresso + ETA para o solver numerico (sem dependencias).

O calculo e dividido em fases ponderadas (matriz, camadas, solve, superficie).
Cada fase reporta uma fracao interna 0..1; a fracao global e a soma dos pesos das
fases ja concluidas mais o peso da fase atual vezes sua fracao. O ETA e
auto-corretivo: eta = decorrido * (1 - fracao) / fracao.

Uso (parametro `progresso` de EstudoNumerico.resolver / resolver_rg):
  - None / False      -> _NullProgresso (no-op, custo zero).
  - callable          -> recebe um dict de evento por passo (uso da GUI).
  - True / "barra" / "cli" -> barra de texto em stderr (so desenha em terminal).

Evento entregue ao callable:
  {"fase", "fracao", "fracao_fase", "decorrido", "eta", "concluido"}
"""

import sys
import time


def _fmt_eta(seg):
    """Formata segundos restantes como '12s' / '1m03s' / '--' (None/invalido)."""
    if seg is None or seg != seg or seg < 0:        # None ou NaN ou negativo
        return "--"
    seg = int(round(seg))
    if seg < 60:
        return f"{seg}s"
    return f"{seg // 60}m{seg % 60:02d}s"


class _NullProgresso:
    """Sem operacao: caminho progresso=None. Todos os metodos custam zero."""

    def iniciar(self):
        pass

    def fase(self, nome):
        pass

    def passo(self, frac_fase):
        pass

    def concluir(self):
        pass


class Progresso:
    """Agrega fases ponderadas e despacha um evento por passo a um sink."""

    def __init__(self, emitir, pesos):
        # descarta fases de peso <= 0 (ex.: 'camadas' em solo uniforme)
        pesos = {k: float(v) for k, v in pesos.items() if v and v > 0}
        total = sum(pesos.values()) or 1.0
        self._pesos = {k: v / total for k, v in pesos.items()}
        self._emitir = emitir
        self._base = 0.0            # peso acumulado das fases ja concluidas
        self._fase = None
        self._peso_atual = 0.0
        self._t0 = None
        self._ultimo = 0.0          # ultima fracao emitida (garante monotonia)

    def iniciar(self):
        self._t0 = time.perf_counter()

    def fase(self, nome):
        if self._fase is not None:                  # fecha a fase anterior
            self._base += self._peso_atual
        self._fase = nome
        self._peso_atual = self._pesos.get(nome, 0.0)

    def passo(self, frac_fase):
        if self._t0 is None:
            self.iniciar()
        frac_fase = 0.0 if frac_fase < 0 else (1.0 if frac_fase > 1 else frac_fase)
        frac = self._base + self._peso_atual * frac_fase
        if frac < self._ultimo:                     # nunca retrocede
            frac = self._ultimo
        self._ultimo = frac
        self._despachar(frac, frac_fase, concluido=False)

    def concluir(self):
        if self._t0 is None:
            self.iniciar()
        self._ultimo = 1.0
        self._despachar(1.0, 1.0, concluido=True)

    def _despachar(self, frac, frac_fase, concluido):
        decorrido = time.perf_counter() - self._t0
        if concluido:
            eta = 0.0
        elif frac <= 0.0:
            eta = None
        else:
            eta = decorrido * (1.0 - frac) / frac
        self._emitir({
            "fase": self._fase,
            "fracao": frac,
            "fracao_fase": frac_fase,
            "decorrido": decorrido,
            "eta": eta,
            "concluido": concluido,
        })


class _BarraTerminal:
    """Sink que desenha uma barra de texto em stderr (carriage-return).

    So desenha quando o fluxo de saida e um terminal (isatty), de modo que pipes
    e redirecionamentos ficam limpos automaticamente. `forcar=True` ignora o teste
    de terminal (util em teste).
    """

    def __init__(self, largura=28, fluxo=None, forcar=False):
        self._largura = largura
        self._fluxo = fluxo if fluxo is not None else sys.stderr
        isatty = getattr(self._fluxo, "isatty", None)
        self._ativo = bool(forcar or (isatty and isatty()))

    def __call__(self, ev):
        if not self._ativo:
            return
        frac = ev["fracao"]
        cheio = int(round(frac * self._largura))
        barra = "#" * cheio + "-" * (self._largura - cheio)
        linha = (f"\r[{barra}] {frac * 100:3.0f}% | "
                 f"{(ev['fase'] or ''):<10} | ETA {_fmt_eta(ev['eta'])}")
        self._fluxo.write(linha)
        if ev["concluido"]:
            self._fluxo.write("\n")
        self._fluxo.flush()


def criar_progresso(progresso, pesos):
    """Fabrica um reporter a partir do parametro `progresso` do usuario.

    None/False -> _NullProgresso; callable -> Progresso(callable);
    True/'barra'/'cli' -> Progresso(_BarraTerminal).
    """
    if progresso is None or progresso is False:
        return _NullProgresso()
    if callable(progresso):
        return Progresso(progresso, pesos)
    if progresso is True or progresso in ("barra", "cli"):
        return Progresso(_BarraTerminal(), pesos)
    raise TypeError("progresso deve ser None, callable, True ou 'barra'/'cli'")
