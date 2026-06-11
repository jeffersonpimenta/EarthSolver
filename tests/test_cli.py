"""Testes de robustez da CLI: entrada invalida vira 'erro: ...' limpo (rc 1)."""

import json

from earthsolver.cli import main

_MALHA = {"area": 100, "Lc": 60, "comprimento_x": 10, "comprimento_y": 10,
          "espac_D": 5, "prof_h": 0.5, "d": 0.01}


def test_analisar_sem_ig_da_erro_limpo(tmp_path, capsys):
    proj = tmp_path / "proj.json"
    proj.write_text(json.dumps({"solo": {"rho": [100.0]},
                                "malha": _MALHA,
                                "falta": {"t": 0.5}}))
    rc = main(["analisar", str(proj)])
    assert rc == 1
    err = capsys.readouterr().err
    assert err.startswith("erro:")
    assert "Ig" in err


def test_estratificar_csv_linha_curta_da_erro_limpo(tmp_path, capsys):
    csv = tmp_path / "sond.csv"
    csv.write_text("spacing,resistance\n1,10\n2\n")
    rc = main(["estratificar", "--entrada", str(csv)])
    assert rc == 1
    assert capsys.readouterr().err.startswith("erro:")
