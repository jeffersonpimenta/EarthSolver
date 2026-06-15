"""Graficos do earthsolver: vista da malha, mapa de potencial de superficie e
mapas de seguranca (tensoes de toque/passo, margem, corrente por segmento).

Usa a API orientada a objeto do matplotlib (`matplotlib.figure.Figure`), sem
`pyplot` nem backend interativo: `savefig` renderiza em PNG via Agg. matplotlib
e importado de forma preguicosa (so quando se plota).
"""

import numpy as np


def _nova_ax():
    from matplotlib.figure import Figure
    fig = Figure(figsize=(6.5, 5.5))
    return fig.add_subplot(111)


def _nova_ax_3d():
    import mpl_toolkits.mplot3d  # noqa: F401  (registra a projecao "3d")
    from matplotlib.figure import Figure
    fig = Figure(figsize=(6.5, 5.5))
    return fig.add_subplot(111, projection="3d")


def _e_haste(c):
    return c.p1[0] == c.p2[0] and c.p1[1] == c.p2[1]


def plot_malha(eletrodo, ax=None):
    """Desenha a vista em planta: condutores como linhas, hastes como marcadores."""
    if ax is None:
        ax = _nova_ax()
    primeiro_cond = primeiro_haste = True
    for c in eletrodo.condutores:
        if _e_haste(c):
            ax.plot([c.p1[0]], [c.p1[1]], marker="v", color="C3", linestyle="",
                    label="haste" if primeiro_haste else None)
            primeiro_haste = False
        else:
            ax.plot([c.p1[0], c.p2[0]], [c.p1[1], c.p2[1]], "-", color="C0",
                    label="condutor" if primeiro_cond else None)
            primeiro_cond = False
    ax.set_aspect("equal")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_title("Malha de aterramento")
    if not (primeiro_cond and primeiro_haste):
        ax.legend(loc="best", fontsize="small")
    return ax


def plot_potencial(raster, eletrodo=None, ax=None):
    """Mapa de potencial de superficie (contourf) a partir de raster (X, Y, Phi).

    Sobrepoe o contorno da malha se `eletrodo` for dado.
    """
    X, Y, Phi = raster
    if ax is None:
        ax = _nova_ax()
    cs = ax.contourf(X, Y, Phi, levels=20)
    ax.figure.colorbar(cs, ax=ax, label="Potencial (V)")
    if eletrodo is not None:
        for c in eletrodo.condutores:
            if _e_haste(c):
                ax.plot([c.p1[0]], [c.p1[1]], marker="v", color="k",
                        markersize=4, linestyle="")
            else:
                ax.plot([c.p1[0], c.p2[0]], [c.p1[1], c.p2[1]], "-",
                        color="k", linewidth=0.6, alpha=0.6)
    ax.set_aspect("equal")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_title("Potencial de superficie")
    return ax


def plot_malha_3d(eletrodo, ax=None):
    """Perspectiva 3D da malha: condutores enterrados e hastes verticais.

    z = profundidade (m); o eixo z e invertido para a profundidade crescer
    para baixo, com a superficie (z=0) no topo.
    """
    if ax is None:
        ax = _nova_ax_3d()
    primeiro_cond = primeiro_haste = True
    for c in eletrodo.condutores:
        xs, ys, zs = ([c.p1[0], c.p2[0]], [c.p1[1], c.p2[1]], [c.p1[2], c.p2[2]])
        if _e_haste(c):
            ax.plot(xs, ys, zs, marker="v", color="C3",
                    label="haste" if primeiro_haste else None)
            primeiro_haste = False
        else:
            ax.plot(xs, ys, zs, "-", color="C0",
                    label="condutor" if primeiro_cond else None)
            primeiro_cond = False
    ax.invert_zaxis()
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_zlabel("profundidade (m)")
    ax.set_title("Malha de aterramento (perspectiva)")
    if not (primeiro_cond and primeiro_haste):
        ax.legend(loc="best", fontsize="small")
    return ax


def plot_potencial_3d(raster, eletrodo=None, ax=None):
    """Elevacao do potencial: superficie 3D Phi(x, y) a partir de raster.

    Desenha a malha projetada no plano da base se `eletrodo` for dado.
    """
    X, Y, Phi = raster
    if ax is None:
        ax = _nova_ax_3d()
    surf = ax.plot_surface(X, Y, Phi, cmap="viridis", linewidth=0,
                           antialiased=True)
    ax.figure.colorbar(surf, ax=ax, label="Potencial (V)", shrink=0.6)
    if eletrodo is not None:
        zbase = float(Phi.min())
        for c in eletrodo.condutores:
            if _e_haste(c):
                ax.plot([c.p1[0]], [c.p1[1]], [zbase], marker="v", color="k",
                        markersize=4, linestyle="")
            else:
                ax.plot([c.p1[0], c.p2[0]], [c.p1[1], c.p2[1]], [zbase, zbase],
                        "-", color="k", linewidth=0.6, alpha=0.6)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_zlabel("Potencial (V)")
    ax.set_title("Potencial de superficie (elevacao)")
    return ax


# --------------------------------------------------------- mapas de seguranca

def _overlay_malha(ax, eletrodo, cor="k"):
    """Sobrepoe o contorno da malha (condutores + hastes) num Axes 2D."""
    if eletrodo is None:
        return
    for c in eletrodo.condutores:
        if _e_haste(c):
            ax.plot([c.p1[0]], [c.p1[1]], marker="v", color=cor,
                    markersize=4, linestyle="")
        else:
            ax.plot([c.p1[0], c.p2[0]], [c.p1[1], c.p2[1]], "-",
                    color=cor, linewidth=0.6, alpha=0.6)


def _bbox_eletrodo(eletrodo):
    """Caixa envolvente (xmin, xmax, ymin, ymax) dos condutores em planta."""
    xs, ys = [], []
    for c in eletrodo.condutores:
        xs += [c.p1[0], c.p2[0]]
        ys += [c.p1[1], c.p2[1]]
    return min(xs), max(xs), min(ys), max(ys)


def _hachura_excedencia(ax, X, Y, campo, limite, bbox=None):
    """Realca onde campo > limite: hachura a regiao e tracа a curva-limite.

    Se `bbox` (xmin, xmax, ymin, ymax) for dado, restringe a hachura a essa
    caixa (a tensao de toque so faz sentido onde ha metal aterrado a tocar).
    """
    campo = np.asarray(campo, dtype=float)
    aval = campo
    if bbox is not None:
        xmin, xmax, ymin, ymax = bbox
        dentro = (X >= xmin) & (X <= xmax) & (Y >= ymin) & (Y <= ymax)
        aval = np.where(dentro, campo, campo.min())     # fora nao excede
    cmax = float(np.nanmax(aval))
    if cmax <= limite:
        return                                           # nada excede o limite
    ax.contourf(X, Y, aval, levels=[limite, cmax], colors="none",
                hatches=["xxx"])
    ax.contour(X, Y, aval, levels=[limite], colors="red", linewidths=1.5)


def _mapa_campo(raster, limite, eletrodo, ax, cmap, rotulo, titulo,
                restringe_bbox):
    """Renderiza um campo escalar (toque ou passo) com realce de excedencia."""
    X, Y, campo = raster
    if ax is None:
        ax = _nova_ax()
    cs = ax.contourf(X, Y, campo, levels=20, cmap=cmap)
    ax.figure.colorbar(cs, ax=ax, label=rotulo)
    if limite is not None:
        bbox = _bbox_eletrodo(eletrodo) \
            if (restringe_bbox and eletrodo is not None) else None
        _hachura_excedencia(ax, X, Y, campo, limite, bbox)
        ax.contour(X, Y, campo, levels=[limite], colors="red",
                   linewidths=1.0, linestyles="--")
    _overlay_malha(ax, eletrodo)
    ax.set_aspect("equal")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_title(titulo)
    return ax


def plot_tensao_toque(raster_toque, limite=None, eletrodo=None, ax=None):
    """Mapa da tensao de toque GPR - Phi(x, y).

    `limite` = tensao de toque toleravel (E_toque): tracа a curva-limite e
    hachura onde excede (restrito a projecao da malha).
    """
    return _mapa_campo(raster_toque, limite, eletrodo, ax, cmap="YlOrRd",
                       rotulo="Tensao de toque (V)",
                       titulo="Tensao de toque", restringe_bbox=True)


def plot_tensao_passo(raster_passo, limite=None, eletrodo=None, ax=None):
    """Mapa da tensao de passo (diferenca de potencial a 1 m, pior direcao).

    `limite` = tensao de passo toleravel (E_passo): a hachura de excedencia
    cobre toda a area (uma pessoa pode pisar em qualquer ponto).
    """
    return _mapa_campo(raster_passo, limite, eletrodo, ax, cmap="PuBuGn",
                       rotulo="Tensao de passo (V)",
                       titulo="Tensao de passo", restringe_bbox=False)


def plot_margem(raster_toque, raster_passo, E_toque, E_passo, eletrodo=None,
                ax=None):
    """Mapa de utilizacao: max(toque/E_toque, passo/E_passo) em %.

    >100 % = reprovado naquele ponto. A contribuicao de toque e restrita a
    projecao da malha; a de passo vale para toda a area.
    """
    X, Y, toque = raster_toque
    _, _, passo = raster_passo
    if ax is None:
        ax = _nova_ax()
    ut = np.asarray(toque, dtype=float) / E_toque
    if eletrodo is not None:
        xmin, xmax, ymin, ymax = _bbox_eletrodo(eletrodo)
        dentro = (X >= xmin) & (X <= xmax) & (Y >= ymin) & (Y <= ymax)
        ut = np.where(dentro, ut, 0.0)
    util = np.maximum(ut, np.asarray(passo, dtype=float) / E_passo) * 100.0
    cs = ax.contourf(X, Y, util, levels=20, cmap="RdYlGn_r")
    ax.figure.colorbar(cs, ax=ax, label="Utilizacao (%)")
    if float(util.max()) > 100.0:
        ax.contourf(X, Y, util, levels=[100.0, float(util.max())],
                    colors="none", hatches=["xxx"])
        ax.contour(X, Y, util, levels=[100.0], colors="black", linewidths=1.5)
    _overlay_malha(ax, eletrodo)
    ax.set_aspect("equal")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_title("Margem de seguranca (utilizacao)")
    return ax


def plot_perfis(raster, raster_toque, raster_passo, E_toque, E_passo, GPR,
                ax=None):
    """Perfis em corte pelas linhas centrais da malha (apresentacao IEEE-80).

    Dois cortes (em x e em y) pelo centro do raster, cada um com potencial de
    superficie, tensao de toque e tensao de passo, mais as linhas de GPR e dos
    limites toleraveis. Devolve a Figure.
    """
    from matplotlib.figure import Figure
    X, Y, Phi = raster
    _, _, toque = raster_toque
    _, _, passo = raster_passo
    ny, nx = Phi.shape
    i, j = ny // 2, nx // 2
    fig = Figure(figsize=(7.0, 8.0))
    cortes = [
        ("Corte em x (y central)", X[i, :], Phi[i, :], toque[i, :], passo[i, :]),
        ("Corte em y (x central)", Y[:, j], Phi[:, j], toque[:, j], passo[:, j]),
    ]
    for k, (titulo, eixo, ph, tq, pa) in enumerate(cortes):
        axk = fig.add_subplot(2, 1, k + 1)
        axk.plot(eixo, ph, color="C0", label="Potencial (V)")
        axk.plot(eixo, tq, color="C1", label="Toque (V)")
        axk.plot(eixo, pa, color="C2", label="Passo (V)")
        axk.axhline(GPR, color="0.4", linestyle=":", linewidth=0.8, label="GPR")
        axk.axhline(E_toque, color="C1", linestyle="--", linewidth=0.8,
                    label="E_toque")
        axk.axhline(E_passo, color="C2", linestyle="--", linewidth=0.8,
                    label="E_passo")
        axk.set_title(titulo)
        axk.set_xlabel("distancia (m)")
        axk.set_ylabel("tensao (V)")
        if k == 0:
            axk.legend(loc="best", fontsize="x-small", ncol=2)
    fig.tight_layout()
    return fig


def plot_corrente(A, B, I, ax=None):
    """Distribuicao da corrente drenada por segmento (planta).

    A, B (M,3): extremos dos segmentos; I (M,): corrente para o solo. Os
    segmentos sao coloridos por |I| (picos tipicamente nos cantos da malha).
    """
    from matplotlib.collections import LineCollection
    if ax is None:
        ax = _nova_ax()
    A = np.asarray(A, dtype=float)
    B = np.asarray(B, dtype=float)
    mag = np.abs(np.asarray(I, dtype=float))
    segs = [[(a[0], a[1]), (b[0], b[1])] for a, b in zip(A, B)]
    lc = LineCollection(segs, array=mag, cmap="plasma", linewidths=3.0)
    ax.add_collection(lc)
    ax.figure.colorbar(lc, ax=ax, label="Corrente drenada |I| (A)")
    # hastes (segmentos verticais) aparecem como ponto em planta: marca-las
    haste = np.isclose(A[:, 0], B[:, 0]) & np.isclose(A[:, 1], B[:, 1])
    if haste.any():
        ax.plot(A[haste, 0], A[haste, 1], "kv", markersize=4, linestyle="")
    ax.autoscale()
    ax.set_aspect("equal")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_title("Distribuicao de corrente por segmento")
    return ax


def salvar(obj, caminho):
    """Salva um Axes ou Figure em PNG. Devolve o caminho."""
    from matplotlib.figure import Figure
    fig = obj if isinstance(obj, Figure) else obj.figure
    fig.savefig(str(caminho), dpi=120, bbox_inches="tight")
    return caminho
