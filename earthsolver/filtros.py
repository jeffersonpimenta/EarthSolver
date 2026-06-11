"""Nucleo numerico do modelo direto de resistividade aparente.

O modelo direto (resistividade aparente de um solo estratificado em N camadas,
para o arranjo de Wenner) e calculado por uma transformada de Hankel. Para manter
o projeto com dependencia unica de numpy (fiel ao EletroSolver), a funcao de
Bessel J0 e avaliada por aproximacao polinomial de Abramowitz & Stegun (cap. 9,
erro < ~1e-7) e a integral e feita por reformulacao numericamente estavel:

    rho_a(a) = rho_1 + 2*a * integral_0^L [T(lmbda) - rho_1] *
                                          [J0(lmbda*a) - J0(2*lmbda*a)] d(lmbda)

onde T(lmbda) e o kernel (transformada de resistividade) obtido pela recursao de
Pekeris. Como [T - rho_1] -> 0 exponencialmente quando lmbda*h_1 cresce, o
integrando tem envelope decrescente e a integral converge sem cauda oscilatoria
lenta (ao contrario da forma com T(lmbda) puro).

Referencias: Abramowitz & Stegun (1972), 9.4.1 e 9.4.3; Koefoed (1979).
"""

import numpy as np

# np.trapz foi renomeado para np.trapezoid no numpy 2.0.
_trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))

# Profundidade adimensional ate a qual [T - rho_1] e desprezivel (e^-24 ~ 4e-11).
_LIM_INTEGRACAO = 12.0
# Numero de pontos da regra do trapezio no intervalo de integracao.
_N_PONTOS = 4001


def j0(x):
    """Funcao de Bessel de primeira especie, ordem 0 (Abramowitz & Stegun)."""
    x = np.asarray(x, dtype=float)
    ax = np.abs(x)
    saida = np.empty_like(ax)

    # Regiao |x| <= 3: serie em t = (x/3)^2.
    pequeno = ax <= 3.0
    t = (x[pequeno] / 3.0) ** 2
    saida[pequeno] = (
        1.0
        - 2.2499997 * t
        + 1.2656208 * t ** 2
        - 0.3163866 * t ** 3
        + 0.0444479 * t ** 4
        - 0.0039444 * t ** 5
        + 0.0002100 * t ** 6
    )

    # Regiao |x| >= 3: forma assintotica com amplitude f0 e fase theta0.
    grande = ~pequeno
    xg = ax[grande]
    z = 3.0 / xg
    f0 = (
        0.79788456
        - 0.00000077 * z
        - 0.00552740 * z ** 2
        - 0.00009512 * z ** 3
        + 0.00137237 * z ** 4
        - 0.00072805 * z ** 5
        + 0.00014476 * z ** 6
    )
    theta0 = (
        xg
        - 0.78539816
        - 0.04166397 * z
        - 0.00003954 * z ** 2
        + 0.00262573 * z ** 3
        - 0.00054125 * z ** 4
        - 0.00029333 * z ** 5
        + 0.00013558 * z ** 6
    )
    saida[grande] = f0 * np.cos(theta0) / np.sqrt(xg)
    return saida


def kernel(lmbda, rho, espessura):
    """Kernel T(lmbda) pela recursao de Pekeris (de baixo para cima).

    rho: resistividades das camadas (rho_1..rho_N).
    espessura: espessuras das N-1 camadas superiores (a ultima e infinita).
    lmbda: array de numeros de onda (1/m).
    """
    lmbda = np.asarray(lmbda, dtype=float)
    T = np.full_like(lmbda, rho[-1], dtype=float)
    for i in range(len(rho) - 2, -1, -1):
        th = np.tanh(lmbda * espessura[i])
        T = (T + rho[i] * th) / (1.0 + (T / rho[i]) * th)
    return T


def resistividade_aparente(rho, espessura, a):
    """Resistividade aparente de Wenner para espacamento(s) `a` (modelo direto).

    Implementa rho_a(a) = rho_1 + 2*a * integral [T - rho_1]*[J0(la)-J0(2la)] dl
    pela regra do trapezio. Aceita `a` escalar ou array.
    """
    rho = np.asarray(rho, dtype=float)
    espessura = np.asarray(espessura, dtype=float)
    rho1 = rho[0]
    escalar = np.isscalar(a) or np.ndim(a) == 0
    avals = np.atleast_1d(np.asarray(a, dtype=float))

    # Limite de integracao governado pela camada mais rasa (decaimento e^-2*l*h1).
    h1 = espessura[0] if espessura.size else 1.0
    limite = _LIM_INTEGRACAO / h1
    lmbda = np.linspace(0.0, limite, _N_PONTOS)
    T_menos_rho1 = kernel(lmbda, rho, espessura) - rho1

    out = np.empty(avals.shape, dtype=float)
    for idx, ai in enumerate(avals):
        g = j0(lmbda * ai) - j0(2.0 * lmbda * ai)
        integrando = T_menos_rho1 * g
        integral = _trapz(integrando, lmbda)
        out[idx] = rho1 + 2.0 * ai * integral

    return float(out[0]) if escalar else out
