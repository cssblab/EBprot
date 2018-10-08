# PluginEBprot

**PluginEBprot** is an implementation of EBprot developed for the [Perseus platform](http://www.coxdocs.org/doku.php?id=perseus:start), providing a GUI for Windows users. 

There are two modules available: 
1. `EBprotV2` runs the differential expression analysis on the given data
2. `EBprotV2.GrpComparisons` contains an extra sub-module that creates ratio data for groups comparisons from multiple experiments and then runs the differential expression analysis like `EBprotV2`

## Full Documentation

See the [Wiki](../../../wiki) for full documentation.

See [Samples](Samples) for sample files and parameters.

A step by step tutorial is available [here](../../../wiki/5.-Getting-Started-with-Perseus-EBprot).
## Requirements

64 bit Windows with .NET Framework 4.5 or higher (See [Perseus Requirements](http://www.coxdocs.org/doku.php?id=perseus:common:download_and_installation))

Latest version of Perseus (See [Perseus Download and Installation Guide](http://www.coxdocs.org/doku.php?id=perseus:common:download_and_installation#download))

## Installing

Upon publication of the EBProt2 paper, the plugin will be available via the Perseus website.

Currently, the plugin is only distributed through github.

To install through github:

* [Download Perseus-PluginEBprot ZIP](../../../releases/download/v1.1.0/pluginEBprot.zip)
* Unzip `pluginEBprot.zip`
* Locate the directory of `Perseus.exe`, which contains `bin` folder
* Copy/Cut `PluginEBprot.dll` file and `EBprotInstallations` folder from `pluginEBprot`
* Paste into the `bin` folder

## Bugs and Feedback

For bugs, questions and discussions please use the [GitHub Issues](../../../issues).
