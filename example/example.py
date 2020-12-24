import pandas as pd
from mppsim import HS_rat_simulation, HS_rat_simulation_random, plot_haplotypes, animate_generations

# Run an HS rat simulation with 30 pairs of rats, circular pair mating, over 60 generations:
sims = HS_rat_simulation(n_pairs=20, n_simulations=1, n_generations=40, chromosomes=[1, 12])
sims.hap_table(generations="all").to_csv("example.txt.gz", sep="\t", index=False)

d = pd.read_csv("example.txt.gz", sep="\t")
d = d[d["chrom"] == 12]

# Plot Chromosome 1 in the population at generation 40:
fig = plot_haplotypes(d[d["generation"] == 40])
fig.savefig("example.png", dpi=300)

# Animate Chromosome 1 in the population across generations:
ani = animate_generations(d, summary=None)
ani.save("example.mp4", dpi=150)
# Plot chromosomes plus haplotype entropy per position:
ani = animate_generations(d, summary="positions")
ani.save("example_positions.mp4", dpi=150)
# Plot chromosomes plus entropy distribution per generation:
ani = animate_generations(d, summary="chrono")
ani.save("example_chrono.mp4", dpi=150)
# Plot chromosomes plus all summary plots:
ani = animate_generations(d, summary="all")
ani.save("example_all.mp4", dpi=150)
