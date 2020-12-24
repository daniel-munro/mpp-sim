.. mppsim documentation master file, created by
   sphinx-quickstart on Thu Jul 30 19:36:03 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

mpp-sim
=======

This package provides tools for genetic simulation of multi-parent populations, i.e. the breeding of diploid organisms over generations. It contains presets for simulating existing model organism populations, though currently only includes HS rats. It contains functions for animated plotting, and will eventually provide tools for comparing real data to simulations in terms of genetic diversity.

Running a simulation
--------------------

A simulation is run using preset specifications for an MPP along with parameters like population size and number of generations. The random assortment and recombination locations of chromosomes are recorded for each generation. Recombinations are determined from the genetic maps for the orgnism's chromosomes.

The basic simulation output is a table giving the spans of founder haplotypes across each chromosome in each individual, in each generation, in each simulation.

Visualization tools
-------------------

Founder haplotypes across chromosomes can be plotted and animated across generations. Diversity (Shannon entropy) per locus can be calculated for a sample of loci for comparison across generations, between breeding methods or parameters, or against diversity calculated from real data.

.. automodule:: mppsim.core
   :members:

.. automodule:: mppsim.plot
   :members:

.. automodule:: mppsim.stats
   :members:

.. toctree::
   :maxdepth: 2
   :caption: Contents:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
