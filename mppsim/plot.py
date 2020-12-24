import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Rectangle
from matplotlib.figure import Figure
from matplotlib.axes import Axes
import matplotlib as mpl
from typing import Iterator
from progress.bar import IncrementalBar

from .stats import entropy_along_chromosome, entropy_summary, region_length_summary


# def region_rectangle(
#     start: float, end: float, y_start: float, height: float, color: tuple
# ) -> Rectangle:
#     width = end - start
#     return Rectangle((start, y_start), width, height, color=color, linewidth=0)


def plot_chromosome(
    table: pd.DataFrame, y_start: float, height: float, colors: dict
) -> Iterator[Rectangle]:
    for i in range(table.shape[0]):
        start = table["start"].iat[i] / 1e6
        end = (table["end"].iat[i] + 1) / 1e6
        yield Rectangle(
            xy=(start, y_start),
            width=end - start,
            height=height,
            color=colors[table["founder"].iat[i]],
            linewidth=0,
        )
        # yield region_rectangle(
        #     start=table["start"].iat[i],
        #     end=table["end"].iat[i],
        #     y_start=y_start,
        #     height=height,
        #     color=color,
        # )


def plot_chromosome_pair(
    table: pd.DataFrame, index: int, gaps: bool, colors: dict
) -> Iterator[Rectangle]:
    d_m = table[table["copy"] == "m"]
    d_p = table[table["copy"] == "p"]
    y_start_m = index + 0.1 if gaps else index
    y_start_p = index + 0.53 if gaps else index + 0.5
    height = 0.37 if gaps else 0.5
    for rect in plot_chromosome(d_m, y_start_m, height, colors):
        yield rect
    for rect in plot_chromosome(d_p, y_start_p, height, colors):
        yield rect


def plot_chromosomes(
    table: pd.DataFrame, gaps: bool, colors: dict
) -> Iterator[Rectangle]:
    # pd.unique preserves order, np.unique does not:
    for i, individual in enumerate(table["individual"].unique()):
        d = table[table["individual"] == individual]
        for rect in plot_chromosome_pair(d, i, gaps, colors):
            yield rect


def plot_haplotypes(table: pd.DataFrame, entropy_panel: bool = True) -> Figure:
    """Plot one chromosome across a population.

    Args:
        table: A table of simulation results containing one chromosome type and one generation.
        entropy_panel: Whether to include a panel for entropy along the chromosome.

    Returns:
        A matplotlib figure.

    """
    for col in ["simulation", "generation", "chromosome"]:
        if col in table:
            assert len(table[col].unique()) == 1
    founders = list("ABCDEFGH")
    colors = {f: plt.get_cmap("Set1")(i) for i, f in enumerate(founders)}
    if entropy_panel:
        fig, (ax, ax2) = plt.subplots(nrows=2, gridspec_kw=dict(height_ratios=[6, 1]))
        # ax = fig.add_subplot()
        entropy = entropy_along_chromosome(table, 500)
        entropy["pos"] = entropy["pos"] / 1e6
        ax2.plot(entropy["pos"], entropy["entropy"])
        ax2.set_ylim((0, np.log2(len(founders))))
        ax2.set_xlabel("Chromosome position (Mb)")
        ax2.set_ylabel("Entropy")
    else:
        fig, ax = plt.subplots()
    # table["start"] = table["start"] / 1e6
    # table["end"] = (table["end"] + 1) / 1e6
    for rect in plot_chromosomes(table, False, colors):
        ax.add_patch(rect)
    ax.set_xlim((1, (table["end"].max() + 1) / 1e6))
    ax.set_ylim((len(table["individual"].unique()), 0))
    if entropy_panel:
        ax2.set_xlim((1, (table["end"].max() + 1) / 1e6))
    else:
        ax.set_xlabel("Chromosome position (Mb)")
    ax.set_ylabel("Individuals")
    return fig


def animate_generation(table: pd.DataFrame, ax: Axes, summary: str):
    founders = list("ABCDEFGH")
    colors = {f: plt.get_cmap("Set1")(i) for i, f in enumerate(founders)}
    # table = table.copy()
    # table["start"] = table["start"] / 1e6
    # table["end"] = (table["end"] + 1) / 1e6
    ax.clear()
    for rect in plot_chromosomes(table, False, colors):
        ax.add_patch(rect)
    # ax.set_yticks([x for x in range(1, len(table["individual"].unique())) if x % 10 == 0])
    ax.set_yticks(
        np.linspace(0, len(table["individual"].unique()), num=6, dtype=int)[1:]
    )
    ax.set_xlim((1, (table["end"].max() + 1) / 1e6))
    ax.set_ylim((len(table["individual"].unique()), 0))
    ax.set_xlabel("Chromosome position (Mb)")
    ax.set_ylabel("Individuals")


def summary_positions(table: pd.DataFrame, ax: Axes):
    founders = list("ABCDEFGH")
    entropy = entropy_along_chromosome(table, 500)
    entropy["pos"] = entropy["pos"] / 1e6
    ax.clear()
    ax.plot(entropy["pos"], entropy["entropy"], color="black")
    ax.set_xlim((1, (table["end"].max() + 1) / 1e6))
    ax.set_yticks(range(int(np.log2(len(founders))) + 2))
    ax.grid(axis="y")
    ax.set_ylim((0, np.log2(len(founders))))
    ax.set_xlabel("Chromosome position (Mb)")
    ax.set_ylabel("Entropy")


def summary_chrono_entropy(stats: pd.DataFrame, ax: Axes, generation: int):
    founders = list("ABCDEFGH")
    last_gen = stats["generation"].max()
    stats = stats[stats["generation"] <= generation]
    ax.clear()
    ax.fill_between(
        x=stats["generation"],
        y1=stats["lower"],
        y2=stats["upper"],
        alpha=0.3,
        color="#555555",
        linewidth=0,
    )
    ax.plot(stats["generation"], stats["median"], color="#555555")
    ax.set_xlim((stats["generation"][0], last_gen))
    ax.set_yticks(range(int(np.log2(len(founders))) + 2))
    ax.grid(axis="y")
    ax.set_ylim((0, np.log2(len(founders))))
    ax.set_xlabel("Generation")
    ax.set_ylabel("Entropy")


def summary_chrono_lengths(stats: pd.DataFrame, ax: Axes, generation: int):
    last_gen = stats["generation"].max()
    min_length = stats["lower"].min()
    max_length = stats["upper"].max()
    stats = stats[stats["generation"] <= generation]
    ax.clear()
    ax.fill_between(
        x=stats["generation"],
        y1=stats["lower"],
        y2=stats["upper"],
        alpha=0.3,
        color="#555555",
        linewidth=0,
    )
    ax.plot(stats["generation"], stats["median"], color="#555555")
    ax.set_xlim((stats["generation"][0], last_gen))
    ax.set_yscale("log", base=10)
    # ax.set_yticks([1e-2, 2e-2, 5e-2, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 1e2, 2e2, 5e2, 1e3])
    ax.set_yticks([1e-3, 1e-2, 0.1, 1, 10, 100, 1000])
    ax.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    ax.grid(axis="y")
    ax.set_ylim((min_length, max_length))
    ax.set_xlabel("Generation")
    ax.set_ylabel("Seg. len. (Mb)")


def animate_generations(
    table: pd.DataFrame, summary: str = "positions", fps: int = 5
) -> FuncAnimation:
    """Create an animated plot of a chromosome across a population over generations.

    Args:
        table: A table of simulation results containing one chromosome type.
        summary: "positions" to show entropy across the chromosome, "chrono" to show entropy distribution and segment length distribution per generation, "all" for all of these, and None for no summary panels.

    Returns:
        A matplotlib.FuncAnimation object.

    """
    def init():
        pass

    def update(g):
        table_g = table[table["generation"] == g]
        animate_generation(table_g, ax, summary)
        if summary in {"positions", "all"}:
            summary_positions(table_g, ax_pos)
        if summary in {"chrono", "all"}:
            summary_chrono_entropy(stats_entropy, ax_ent, g)
            summary_chrono_lengths(stats_lengths, ax_len, g)
        ax.text(
            0.5,
            1.05,
            "Generation {}".format(g),
            size=plt.rcParams["axes.titlesize"],
            ha="center",
            transform=ax.transAxes,
        )
        bar.next()
        if g == gens[-1]:
            bar.finish()

    for col in ["simulation", "chromosome"]:
        if col in table:
            assert len(table[col].unique()) == 1
    if summary == "positions":
        fig, (ax, ax_pos) = plt.subplots(
            nrows=2, gridspec_kw=dict(height_ratios=[6, 1], hspace=0.3), figsize=(8, 6)
        )
    elif summary == "chrono":
        fig, axs = plt.subplots(
            nrows=2,
            ncols=2,
            gridspec_kw=dict(height_ratios=[6, 1], hspace=0.3, wspace=0.4),
            figsize=(8, 6),
        )
        axs[0, 0].remove()
        axs[0, 1].remove()
        ax = fig.add_subplot(axs[0, 0].get_gridspec()[0, :])
        ax_ent = axs[1, 0]
        ax_len = axs[1, 1]
        # Store 5%, 50%, and 95% entropy/length quantiles per generation:
        stats_entropy = entropy_summary(table, 200)
        stats_lengths = region_length_summary(table)
    elif summary == "all":
        fig, axs = plt.subplots(
            nrows=3,
            ncols=2,
            gridspec_kw=dict(height_ratios=[5, 1, 1], hspace=0.6, wspace=0.4),
            figsize=(8, 6),
        )
        for r, c in [(0, 0), (0, 1), (1, 0), (1, 1)]:
            axs[r, c].remove()
        ax = fig.add_subplot(axs[0, 0].get_gridspec()[0, :])
        ax_pos = fig.add_subplot(axs[1, 0].get_gridspec()[1, :])
        ax_ent = axs[2, 0]
        ax_len = axs[2, 1]
        # Store 5%, 50%, and 95% entropy/length quantiles per generation:
        stats_entropy = entropy_summary(table, 200)
        stats_lengths = region_length_summary(table)
    else:
        fig, ax = plt.subplots(figsize=(8, 6))
    gens = np.sort(table["generation"].unique())
    bar = IncrementalBar("Generations plotted", max=len(gens))
    ani = FuncAnimation(fig, update, frames=gens, init_func=init, interval=1000 // fps)
    return ani

