# CLI EBprot

**CLI EBprot** is a C++ implementation of EBprot offering flexible semi-parametric modeling on top of the [version 1](http://onlinelibrary.wiley.com/doi/10.1002/pmic.201400620/abstract;jsessionid=613BD152847535F4278E0C6ED6ACD036.f02t01) (Koh et al, *Proteomics* 2015) implemented in R.

There are two modules available: 
1. `EBprotV2` is the main module that performs a mixture modeling on the log ratio data and estimates the posterior probability of differential expression for a protein and the corresponding FDR.

2. `EBprot.MakeGrpData` is a sub-module which helps to create ratio data for group comparisons from multiple experiments. This generates a weighted averate ratio for differential expression analysis using EBprotV2. This module should be used before EBprotV2 in datasets where there are multiple samples across different groups. 

## Full Documentation

See the [Wiki](../../../wiki) for full documentation.

See [Example](Example) for sample files.

A step by step tutorial is available [here](../../../wiki/1.-Getting-Started-with-CLI-EBprot).

## Requirements

GCC Compiler

GNU Scientific Library (GSL) which can be downloaded [here](https://www.gnu.org/software/gsl/)

[R programming language](https://www.r-project.org) to visualize the outputs of the program

## Installing

* [Download EBprotV2 ZIP](../../../releases/download/v1.1.0/EBprotV2.zip)
* Type `make` to install all the components of the software

## Bugs and Feedback

For bugs, questions and discussions please use the [GitHub Issues](../../../issues).

