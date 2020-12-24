from bisect import bisect_left, bisect_right
from collections import OrderedDict
import numpy as np
import pandas as pd


def physical_pos(gen_map, gen_pos):
    r = bisect_left(gen_map["cm"], gen_pos)
    if r == len(gen_map["cm"]):
        return gen_map["pos"][r - 1]
    elif gen_map["cm"][r] == gen_pos or r == 0:
        return gen_map["pos"][r]
    else:
        # Interpolate the physical position.
        g_lo = gen_map["cm"][r - 1]
        g_hi = gen_map["cm"][r]
        p_lo = gen_map["pos"][r - 1]
        p_hi = gen_map["pos"][r]
        rel = (gen_pos - g_lo) / (g_hi - g_lo)
        return p_lo + round(rel * (p_hi - p_lo))


def load_gen_maps(format_string: str, chromosomes: list) -> OrderedDict:
    gen_maps = OrderedDict()
    for chrom in chromosomes:
        filename = format_string.format(chrom)
        gen_maps[chrom] = pd.read_table(filename, sep=" ", names=["pos", "ratio", "cm"])
    return gen_maps


class Chromosome:
    def __init__(self, name: str, gen_map: OrderedDict):
        self.name = name
        self.gen_map = gen_map
        self.gen_length = gen_map["cm"].iloc[-1]
        self.phys_length = gen_map["pos"].iloc[-1]

    def initialize(self, founder: str):
        self.haps = [[founder], [founder]]
        self.starts = [[1], [1]]

    def crossovers(self) -> list:
        rng = np.random.default_rng()
        n_crossovers = rng.poisson(self.gen_length / 100)
        gen_locs = np.sort(rng.uniform(0, self.gen_length, size=n_crossovers))
        crosses = [physical_pos(self.gen_map, loc) for loc in gen_locs]
        # Remove any "crossover" at base 1:
        return [cross for cross in crosses if cross > 1]

    def meiosis(self) -> tuple:
        crosses = self.crossovers()
        h = np.random.default_rng().choice(2)
        if len(crosses) == 0:
            return (self.starts[h][:], self.haps[h][:])
        # If there are crosses, copy regions while alternating between chromosomes.
        i = bisect_left(self.starts[h], crosses[0])
        starts, haps = self.starts[h][:i], self.haps[h][:i]
        # Add regions between crossovers.
        for c, cross in enumerate(crosses[:-1]):
            h = (h + 1) % 2
            # Add haplotype at crossover if different from previous.
            cross_hap = self.haps[h][bisect_left(self.starts[h], cross) - 1]
            if haps[-1] != cross_hap:
                starts.append(cross)
                haps.append(cross_hap)
            # Then, add haplotypes between this and next cross.
            i = bisect_right(self.starts[h], cross)
            j = bisect_left(self.starts[h], crosses[c + 1])
            starts.extend(self.starts[h][i:j])
            haps.extend(self.haps[h][i:j])
        h = (h + 1) % 2
        # Add haplotype at last crossover if different from previous.
        cross_hap = self.haps[h][bisect_left(self.starts[h], crosses[-1]) - 1]
        if haps[-1] != cross_hap:
            starts.append(crosses[-1])
            haps.append(cross_hap)
        # Then, add haplotypes after last cross.
        i = bisect_right(self.starts[h], crosses[-1])
        starts.extend(self.starts[h][i:])
        haps.extend(self.haps[h][i:])
        return (starts, haps)

    def hap_table(self) -> pd.DataFrame:
        ends1 = [s - 1 for i, s in enumerate(self.starts[0][1:])]
        ends1.append(self.phys_length)
        ends2 = [s - 1 for i, s in enumerate(self.starts[1][1:])]
        ends2.append(self.phys_length)
        haps = pd.DataFrame(
            {
                "chrom": self.name,
                "copy": ["m"] * len(ends1) + ["p"] * len(ends2),
                "start": self.starts[0] + self.starts[1],
                "end": ends1 + ends2,
                "founder": self.haps[0] + self.haps[1],
            }
        )
        return haps


class Individual:
    def __init__(self, gen_maps: OrderedDict):
        self.gen_maps = gen_maps
        self.chromosomes = []
        for chrom in gen_maps.keys():
            chromosome = Chromosome(chrom, gen_maps[chrom])
            self.chromosomes.append(chromosome)

    def initialize(self, founder: str):
        for chromosome in self.chromosomes:
            chromosome.initialize(founder)

    def offspring(self, father: "Individual") -> "Individual":
        child = Individual(self.gen_maps)
        for i, chrom_m in enumerate(self.chromosomes):
            chrom = child.chromosomes[i]
            starts_m, haps_m = chrom_m.meiosis()
            starts_p, haps_p = father.chromosomes[i].meiosis()
            chrom.starts = [starts_m, starts_p]
            chrom.haps = [haps_m, haps_p]
        return child

    def hap_table(self) -> pd.DataFrame:
        return pd.concat([chrom.hap_table() for chrom in self.chromosomes])


class Population:
    def __init__(self, n_pairs: int):
        self.pairs = [{} for i in range(n_pairs)]

    def initialize(self, gen_maps: OrderedDict, founders: list):
        for pair in self.pairs:
            rng = np.random.default_rng()
            founder = rng.choice(founders, 2, replace=False)
            pair["female"] = Individual(gen_maps)
            pair["female"].initialize(founder[0])
            pair["male"] = Individual(gen_maps)
            pair["male"].initialize(founder[1])

    def next_generation(self, mating_scheme: str) -> "Population":
        new_gen = Population(len(self.pairs))
        if mating_scheme == "circular":
            i_father = [len(self.pairs) - 1] + list(range(len(self.pairs) - 1))
        elif mating_scheme == "random":
            # Permute fathers until no pairs are siblings (i.e. find a derangement):
            rng = np.random.default_rng()
            while True:
                i_father = rng.permutation(len(self.pairs))
                if not np.any(i_father == range(len(self.pairs))):
                    break
        for i, pair in enumerate(new_gen.pairs):
            mother = self.pairs[i]["female"]
            father = self.pairs[i_father[i]]["male"]
            pair["female"] = mother.offspring(father)
            pair["male"] = mother.offspring(father)
        return new_gen

    def hap_table(self) -> pd.DataFrame:
        haps = []
        for i, pair in enumerate(self.pairs):
            female = pair["female"].hap_table()
            female["individual"] = 2 * i + 1
            haps.append(female)
            male = pair["male"].hap_table()
            male["individual"] = 2 * i + 2
            haps.append(male)
        haps = pd.concat(haps)
        return haps[["individual", "chrom", "copy", "start", "end", "founder"]]


class Simulation:
    def __init__(self, founders: list, mating_scheme: str, n_pairs: int, gen_maps: OrderedDict):
        pop_0 = Population(n_pairs)
        pop_0.initialize(gen_maps, founders)
        self.generations = [pop_0]
        self.mating_scheme = mating_scheme

    def add_generations(self, n_generations: int = 1):
        for i in range(n_generations):
            # print("Simulating Generation {}...".format(len(self.generations)))
            new_gen = self.generations[-1].next_generation(self.mating_scheme)
            self.generations.append(new_gen)

    def hap_table(self, generations="last") -> pd.DataFrame:
        haps = []
        if generations == "last":
            generations = [len(self.generations) - 1]
        elif generations == "all":
            generations = range(len(self.generations))
        for i in generations:
            gen = self.generations[i].hap_table()
            gen["generation"] = i
            haps.append(gen)
        haps = pd.concat(haps)
        return haps[
            ["generation", "individual", "chrom", "copy", "start", "end", "founder"]
        ]


class Simulations:
    def __init__(self, founders: list, mating_scheme: str, n_pairs: int, gen_maps: OrderedDict, n_simulations: int, n_generations: int):
        self.simulations = []
        for i in range(n_simulations):
            print("Running simulation {}...".format(i + 1))
            simulation = Simulation(founders, mating_scheme, n_pairs, gen_maps)
            simulation.add_generations(n_generations)
            self.simulations.append(simulation)

    def hap_table(self, generations="last") -> pd.DataFrame:
        haps = []
        for i, simulation in enumerate(self.simulations):
            sim = simulation.hap_table(generations)
            sim["simulation"] = i + 1
            haps.append(sim)
        haps = pd.concat(haps)
        return haps[
            [
                "simulation",
                "generation",
                "individual",
                "chrom",
                "copy",
                "start",
                "end",
                "founder",
            ]
        ]


def HS_rat_simulation(n_pairs: int, n_simulations: int, n_generations: int, chromosomes: list = range(1, 21)) -> Simulations:
    """Simulate an HS rat population with circular pair mating.

    Args:
        n_pairs: Number of pairs of rats in the population (i.e. population size is twice this).
        n_simulations: Number of independent simulations to run.
        n_generations: Number of generations to run per simulation (not including the initial population stage).
        chromosomes: List of chromosomes (as integers) to simulate.

    Returns:
        A Simulations object.

    """
    gen_maps = load_gen_maps("~/eqtls/data/genotype/genetic_map/MAP4chr{}.txt", chromosomes)
    founders = list("ABCDEFGH")
    return Simulations(founders, "circular", n_pairs, gen_maps, n_simulations, n_generations)


def HS_rat_simulation_random(n_pairs: int, n_simulations: int, n_generations: int, chromosomes: list = range(1, 21)) -> Simulations:
    """Simulate an HS rat population with random pair mating.

    Args:
        n_pairs: Number of pairs of rats in the population (i.e. population size is twice this).
        n_simulations: Number of independent simulations to run.
        n_generations: Number of generations to run per simulation (not including the initial population stage).
        chromosomes: List of chromosomes (as integers) to simulate.

    Returns:
        A Simulations object.

    """
    gen_maps = load_gen_maps("~/eqtls/data/genotype/genetic_map/MAP4chr{}.txt", chromosomes)
    founders = list("ABCDEFGH")
    return Simulations(founders, "random", n_pairs, gen_maps, n_simulations, n_generations)
