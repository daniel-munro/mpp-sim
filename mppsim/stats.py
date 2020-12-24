from bisect import bisect_right
import numpy as np
import pandas as pd
from progress.bar import IncrementalBar


def haplotype_at_position(table: pd.DataFrame, pos: int) -> str:
    indices = table["start"].searchsorted(pos, "right") - 1
    return table["founder"].iloc[indices]


def entropy(counts: list) -> float:
    p = np.array(counts) / sum(counts)
    return -sum(p[np.nonzero(p)] * np.log2(p[np.nonzero(p)]))


def entropy_along_chromosome(table: pd.DataFrame, n_pos: int) -> pd.DataFrame:
    assert len(table["chrom"].unique()) == 1
    posns = np.linspace(1, table["end"].max(), num=n_pos, dtype=int)
    codes, uniques = table["founder"].factorize()
    table = table.copy()
    table["founder"] = codes
    grouped = table.groupby(["simulation", "generation", "individual", "chrom", "copy"])
    haps = np.empty((len(grouped), n_pos), dtype=int)
    for i, (name, group) in enumerate(grouped):
        haps[i, :] = haplotype_at_position(group, posns)
    entropies = [entropy(np.unique(haps[:, i], return_counts=True)[1]) for i in range(n_pos)]
    return pd.DataFrame({"pos": posns, "entropy": entropies})


def haplotype_region_lengths(table: pd.DataFrame) -> np.ndarray:
    return np.sort(table["end"] - table["start"] + 1) / 1e6


def entropy_summary(table: pd.DataFrame, n_pos: int) -> pd.DataFrame:
    """Store 5%, 50%, and 95% entropy quantiles per generation.

    Args:
        table: A table of simulation results containing one chromosome type.
        n_pos: Number of equally-spaced positions at which to sample entropy.

    Returns:
        A table with columns "generation", "lower", "median", and "upper".
    
    """
    stats = dict(generation=[], lower=[], median=[], upper=[])
    grouped = table.groupby("generation")
    bar = IncrementalBar("Generations summarized", max=len(grouped))
    for g, table_g in grouped:
        entropy = entropy_along_chromosome(table_g, n_pos)["entropy"]
        quantiles = np.quantile(entropy, [0.05, 0.5, 0.95])
        stats["generation"].append(g)
        stats["lower"].append(quantiles[0])
        stats["median"].append(quantiles[1])
        stats["upper"].append(quantiles[2])
        bar.next()
    bar.finish()
    return pd.DataFrame(stats)


def region_length_summary(table: pd.DataFrame) -> pd.DataFrame:
    """Store 5%, 50%, and 95% length quantiles per generation.
    
    Args:
        table: A table of simulation results containing one chromosome type.

    Returns:
        A table with columns "generation", "lower", "median", and "upper".
    
    """
    grouped = table.groupby("generation")
    stats = dict(generation=[], lower=[], median=[], upper=[])
    for g, table_g in grouped:
        lengths = haplotype_region_lengths(table_g)
        quantiles = np.quantile(lengths, [0.05, 0.5, 0.95])
        stats["generation"].append(g)
        stats["lower"].append(quantiles[0])
        stats["median"].append(quantiles[1])
        stats["upper"].append(quantiles[2])
    return pd.DataFrame(stats)
