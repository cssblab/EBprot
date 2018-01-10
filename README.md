# PluginEBprot

Perseus Plugin for EBprot

**EBprot** is a statistical tool for differential expression of proteins in labeling-based proteomics experiments. The tool directly models the distribution of peptide-level ratios and scores each protein based on consistent evidence across the peptides, in contrast to the common practice of assigning p-values to protein-level ratios which loses reproducibility information.

The modules presented here contain an improved extension of EBprot (version 1) offering flexible semi-parametric modeling which is more robust in handling heavily skewed or multi-modal data, like ratio data. 

The principal method (version 1) was published [here](http://onlinelibrary.wiley.com/doi/10.1002/pmic.201400620/abstract;jsessionid=613BD152847535F4278E0C6ED6ACD036.f02t01) (Koh et al, *Proteomics* 2015). 

## Full Documentation

See the [Wiki](../../wiki) for full documentation.

A step by step tutorial is available [here]().

## Requirements

64 bit Windows with .NET Framework 4.5 or higher (See [Perseus Requirements](http://www.coxdocs.org/doku.php?id=perseus:common:download_and_installation))

PERSEUS version 1.6.0.2 (See [Perseus Download and Installation Guide](http://www.coxdocs.org/doku.php?id=perseus:common:download_and_installation#download))

## Installing

Upon publication of the EBProt2 paper, the plugin will be available via the Perseus website.

Currently, the plugin is only distributed through github.

To install through github:

* [Download Perseus-PluginEBprot ZIP](../../releases/download/v0.0.1/pluginEBprot_0.0.1.zip)
* Unzip `pluginEBport_<version>.zip`
* Locate the directory of `Perseus.exe`, which contains `bin` folder
* Copy/Cut `PluginEBprot.dll` file and `EBprotInstallations` folder from `pluginEBprot_<version>.zip\Plugin`
* Paste into the `bin` folder

## Bugs and Feedback

For bugs, questions and discussions please use the [GitHub Issues](../../issues).
