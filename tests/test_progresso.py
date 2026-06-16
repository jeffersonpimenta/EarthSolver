"""Testes da barra de progresso / ETA do solver numerico.

Verifica o hook `progresso=` (eventos monotonicos, fases observadas, conclusao),
a ausencia de regressao no caminho sem progresso e os utilitarios do modulo.
"""

import io

import pytest

from earthsolver.numerico import (Condutor, Eletrodo, EstudoNumerico,
                                  estudo_convergencia)
from earthsolver.progresso import (Progresso, _BarraTerminal, _NullProgresso,
                                   _fmt_eta, criar_progresso)
from earthsolver.solo import ModeloSolo


def _malha(**kw):
    return Eletrodo.malha_retangular(70.0, 70.0, 7.0, 7.0, 0.5, 0.01, **kw)


def _haste(L=3.0):
    return Eletrodo([Condutor((0.0, 0.0, 0.0), (0.0, 0.0, L), 0.01)])


# ----------------------------------------------------------------- integracao
def test_callback_recebe_eventos_monotonicos():
    solo = ModeloSolo([400.0], [])
    eventos = []
    est = EstudoNumerico(solo, _malha(), Ig=1908.0, t=0.5, comp_alvo=7.0)
    est.resolver(progresso=eventos.append)

    assert eventos
    fracs = [e["fracao"] for e in eventos]
    assert fracs == sorted(fracs)                       # monotonica nao-decrescente
    assert all(0.0 <= f <= 1.0 for f in fracs)
    assert all(0.0 <= e["fracao_fase"] <= 1.0 for e in eventos)

    ultimo = eventos[-1]
    assert ultimo["concluido"] is True
    assert ultimo["fracao"] == pytest.approx(1.0)
    assert ultimo["eta"] == 0.0

    fases = {e["fase"] for e in eventos}
    assert {"matriz", "solve", "superficie"} <= fases


def test_camadas_aparece_em_solo_multicamada():
    solo = ModeloSolo([100.0, 400.0], [1.5])
    eventos = []
    est = EstudoNumerico(solo, _haste(), Ig=10.0, t=0.5, comp_alvo=0.5)
    est.resolver(progresso=eventos.append)
    assert "camadas" in {e["fase"] for e in eventos}


def test_progresso_none_nao_altera_resultado():
    solo = ModeloSolo([400.0], [])
    a = EstudoNumerico(solo, _malha(), Ig=1908.0, t=0.5, comp_alvo=7.0).resolver()
    b = EstudoNumerico(solo, _malha(), Ig=1908.0, t=0.5,
                       comp_alvo=7.0).resolver(progresso=lambda e: None)
    for k in ("Rg", "Em", "Es", "GPR"):
        assert a[k] == pytest.approx(b[k])


def test_prog_resetado_apos_resolver():
    solo = ModeloSolo([400.0], [])
    est = EstudoNumerico(solo, _malha(), Ig=1908.0, t=0.5, comp_alvo=7.0)
    est.resolver(progresso=lambda e: None)
    assert isinstance(est._prog, _NullProgresso)        # nao vaza p/ chamadas seguintes


def test_resolver_rg_direto_com_callback():
    solo = ModeloSolo([400.0], [])
    eventos = []
    est = EstudoNumerico(solo, _malha(), Ig=1908.0, t=0.5, comp_alvo=7.0)
    est.resolver_rg(progresso=eventos.append)
    assert eventos and eventos[-1]["concluido"] is True
    assert isinstance(est._prog, _NullProgresso)


def test_convergencia_com_callback():
    solo = ModeloSolo([400.0], [])
    eventos = []
    estudo_convergencia(solo, _malha(), 1908.0, 0.5, [14.0, 7.0],
                        progresso=eventos.append)
    assert eventos
    assert eventos[-1]["concluido"] is True
    assert eventos[-1]["fracao"] == pytest.approx(1.0)


# ----------------------------------------------------------------- unidade
def test_fmt_eta():
    assert _fmt_eta(None) == "--"
    assert _fmt_eta(float("nan")) == "--"
    assert _fmt_eta(-1) == "--"
    assert _fmt_eta(5) == "5s"
    assert _fmt_eta(63) == "1m03s"


def test_criar_progresso_tipos():
    assert isinstance(criar_progresso(None, {}), _NullProgresso)
    assert isinstance(criar_progresso(False, {}), _NullProgresso)
    assert isinstance(criar_progresso(lambda e: None, {"a": 1.0}), Progresso)
    assert isinstance(criar_progresso("barra", {"a": 1.0}), Progresso)
    with pytest.raises(TypeError):
        criar_progresso(123, {})


def test_barra_terminal_desenha_em_tty():
    class _TTY(io.StringIO):
        def isatty(self):
            return True

    buf = _TTY()
    prog = Progresso(_BarraTerminal(fluxo=buf), {"x": 1.0})
    prog.iniciar()
    prog.fase("x")
    prog.passo(0.5)
    prog.concluir()
    saida = buf.getvalue()
    assert "%" in saida and "ETA" in saida
    assert saida.endswith("\n")


def test_barra_terminal_silenciosa_sem_tty():
    buf = io.StringIO()                                 # StringIO.isatty() -> False
    prog = Progresso(_BarraTerminal(fluxo=buf), {"x": 1.0})
    prog.iniciar()
    prog.fase("x")
    prog.passo(0.5)
    prog.concluir()
    assert buf.getvalue() == ""                         # nada quando nao e terminal
