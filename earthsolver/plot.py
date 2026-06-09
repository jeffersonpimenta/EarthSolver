"""Graficos do earthsolver: vista da malha e mapa de potencial de superficie.

Usa a API orientada a objeto do matplotlib (`matplotlib.figure.Figure`), sem
`pyplot` nem backend interativo: `savefig` renderiza em PNG via Agg. matplotlib
e importado de forma preguicosa (so quando se plota).
"""


def _nova_ax():
    from matplotlib.figure import Figure
    fig = Figure(figsize=(6.5, 5.5))
    return fig.add_subplot(111)


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


def salvar(obj, caminho):
    """Salva um Axes ou Figure em PNG. Devolve o caminho."""
    from matplotlib.figure import Figure
    fig = obj if isinstance(obj, Figure) else obj.figure
    fig.savefig(str(caminho), dpi=120, bbox_inches="tight")
    return caminho
