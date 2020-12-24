# mpp-sim: multi-parent population simulator

This package provides tools for genetic simulation of multi-parent populations, i.e. the breeding of diploid organisms over generations. It contains presets for simulating existing model organism populations, though currently only includes HS rats. It contains functions for animated plotting, and will eventually provide tools for comparing real data to simulations in terms of genetic diversity.

## Dependencies

- Python 3.6+
- numpy
- pandas
- progress
- matplotlib

## Installation

Clone the repository or [download and unzip](https://github.com/daniel-munro/mpp-sim/archive/master.zip).  To get Python dependencies and be able to run or import the modules from any directory, install with pip (or pip3 if pip = pip2)::

```
pip install -e mpp-sim/
```

## Running a simulation

See the documentation at docs/_build/html/index.html.

A simulation is run using preset specifications for an MPP along with parameters like population size and number of generations. The random assortment and recombination locations of chromosomes are recorded for each generation. Recombinations are determined from the genetic maps for the orgnism's chromosomes.

The basic simulation output is a table giving the spans of founder haplotypes across each chromosome in each individual, in each generation, in each simulation.

## Visualization tools

Founder haplotypes across chromosomes can be plotted and animated across generations. Diversity (Shannon entropy) per locus can be calculated for a sample of loci for comparison across generations, between breeding methods or parameters, or against diversity calculated from real data.

## TODO

- Add more organisms and breeding schemes.
- Add fine-grained control like F/M ratio and changes in breeding schemes over time.
- Use founder genotypes to output simulated population genotypes.
- Allow easy comparison of genetic diversity across MPPs, breeding schemes, parameters, etc.
- Add feature to compare simulated genetic diversity to real data.
- Add support for sex chromosomes.
