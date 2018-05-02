# PluginEBprot

**EBprot** is a statistical tool for differential expression of proteins in labeling-based proteomics experiments. The tool directly models the distribution of peptide-level ratios and scores each protein based on consistent evidence across the peptides, in contrast to the common practice of assigning p-values to protein-level ratios which loses reproducibility information.

The modules presented here contain an improved extension of EBprot (version 1) offering flexible semi-parametric modeling which is more robust in handling heavily skewed or multi-modal data, like ratio data. 

There are two modules available: 
1. `EBprotV2` runs the differential expression analysis on the given data
2. `EBprotV2.GrpComparisons` contains an extra sub-module that creates ratio data for groups comparisons from multiple experiments and then runs the differential expression analysis like `EBprotV2`

The principal method (version 1) was published [here](http://onlinelibrary.wiley.com/doi/10.1002/pmic.201400620/abstract;jsessionid=613BD152847535F4278E0C6ED6ACD036.f02t01) (Koh et al, *Proteomics* 2015). 

## Full Documentation

See the [Wiki](../../wiki) for full documentation.

A step by step tutorial is available [here](../../wiki/Getting-Started).

## Bugs and Feedback

For bugs, questions and discussions please use the [GitHub Issues](../../issues).

