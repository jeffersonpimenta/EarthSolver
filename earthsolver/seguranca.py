"""Tensoes toleraveis de toque e passo (IEEE Std 80 / NBR 15751).

Formula compartilhada pelo metodo simplificado (malha.py) e pelo solver numerico
(numerico.py): fator de reducao da camada superficial Cs e os limites de toque e
passo para o peso corporeo de referencia (50 ou 70 kg).
"""

import math

# Constante corporea c do IEEE Std 80 (0.116 p/ 50 kg, 0.157 p/ 70 kg).
C_PESO = {50: 0.116, 70: 0.157}


def fator_cs(rho, rho_s, h_s):
    """Fator de reducao da camada superficial Cs (IEEE Std 80).

    rho: resistividade do solo (Ohm.m); rho_s: da camada superficial (brita);
    h_s: espessura da camada superficial (m). Cs = 1 quando rho_s == rho.
    """
    return 1.0 - 0.09 * (1.0 - rho / rho_s) / (2.0 * h_s + 0.09)


def tensoes_toleraveis(rho, rho_s, h_s, t, peso):
    """Tensoes de toque e passo toleraveis (V). Retorna (E_toque, E_passo, Cs).

    rho, rho_s, h_s: solo / camada superficial; t: duracao do choque (s);
    peso: 50 ou 70 (kg).
    """
    if peso not in C_PESO:
        raise ValueError("peso deve ser 50 ou 70 (kg)")
    Cs = fator_cs(rho, rho_s, h_s)
    c = C_PESO[peso]
    E_toque = (1000.0 + 1.5 * Cs * rho_s) * c / math.sqrt(t)
    E_passo = (1000.0 + 6.0 * Cs * rho_s) * c / math.sqrt(t)
    return E_toque, E_passo, Cs
