# mlscluster

[![DOI](https://zenodo.org/badge/518473195.svg)](https://zenodo.org/doi/10.5281/zenodo.10520061)

**[m]ulti[l]evel [s]election [cluster]ing**: tree-based variant fitness clustering algorithm

*mlscluster* is an algorithm that computes three statistics (ratio of descendant 
clade sizes, ratio of persistence times, and logistic growth rates) for sister clades 
in huge time-scaled phylogenies. Clade-defining homoplasies are matched with these 
statistics and results for different cluster thresholds are produced. 
Therefore, this method highlights mutations (or transmission fitness polymorphisms [TFPs]) 
which likely present different selective pressures within and between hosts (multilevel selection). 
This package contains functions to:

 * Compute the multilevel selectives pressures across thresholds;
 * Perform statistical analysis and graphical visualisation for all thresholds to identify lineages 
and genomic regions enriched for TFPs;
 * Quantify the false discovery rate of the algorithm for each threshold;
 * Generate multiple lineage and gene-specific plots of TFPs for the threshold selected by the user.

For a complete description of the statistical methods incorporated into this
package, see our paper:

**Phylogenetic signatures reveal multilevel selection and fitness costs in SARS-CoV-2**. 
Vinicius B. Franceschi & Erik M. Volz. 
[...]

## Installation

In R, install the `devtools` package and run:

```
devtools::install_github('mrc-ide/mlscluster')
```

## Documentation

See the [SARS-CoV-2 computational experiments from our paper](https://github.com/vinibfranc/mlscluster-experiments) for examples of how to use *mlscluster*.

If you need guidance on how to run some functions, see their help pages. For example, within R you can try: 

```
help(mlsclust)
```

## Maintenance

The package currently support SARS-CoV-2 analyses only. Further support for other pathogens, especially HIV, is planned. This involves adding genomic coordinates and specific sanity checks for known sequencing artifacts of each pathogen. Watch this space for future releases or contact us if you are interested in such applications.