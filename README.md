# EBprot

**EBprot** is a statistical tool for differential expression of proteins in labeling-based proteomics experiments. The tool directly models the distribution of peptide-level ratios and scores each protein based on consistent evidence across the peptides, in contrast to the common practice of assigning p-values to protein-level ratios which loses reproducibility information.

The modules presented here contain an improved extension of EBprot (version 1) offering flexible semi-parametric modeling which is more robust in handling heavily skewed or multi-modal data, like ratio data. 

The principal method (version 1) was published [here](http://onlinelibrary.wiley.com/doi/10.1002/pmic.201400620/abstract;jsessionid=613BD152847535F4278E0C6ED6ACD036.f02t01) (Koh et al, *Proteomics* 2015). 

This tool is available as [C++](c%2B%2B) implementation and as a plugin for [Perseus](perseus-plugin) (GUI for Windows users only).

## Full Documentation

See the [Wiki](../../wiki) for full documentation.

See [C++](c%2B%2B) for C++ implementation.

See [Perseus-Plugin](perseus-plugin) for the plugin implementation of Perseus software package.

## Bugs and Feedback

For bugs, questions and discussions please use the [GitHub Issues](../../issues).

## License

Copyright (C) <2018> Hiromi W.L. Koh < hiromi_koh@nuhs.edu.sg >, Yunbin Zhang < yz2236@nyu.edu >, Christine Vogel < cvogel@nyu.edu >, and Hyungwon Choi < hwchoi@nus.edu.sg >, National University of Singapore.

Licensed under the Apache License, Version 2.0 (the "License");

you may not use this file except in compliance with the License.

You may obtain a copy of the License at

[Apache 2.0 license](http://www.apache.org/licenses/LICENSE-2.0)

Unless required by applicable law or agreed to in writing, software

distributed under the License is distributed on an "AS IS" BASIS,

WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.

See the License for the specific language governing permissions and

limitations under the License.


