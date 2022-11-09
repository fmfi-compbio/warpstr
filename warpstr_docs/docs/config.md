# Configuration file

The only required parameter for running WarpSTR is the path to the config file. The file must be in YAML format. It must be populated with elements such as `inputs`, `output` and `reference_path`. An example is provided in `example/config.yaml`.

There are also many advanced parameters that are optional to set. List of all parameters are found in `example/advanced_params.yaml`. To set values for those parameters, just add those parameters to your main config and set them to the desired value. In other case, default values for those parameters are taken.

## Description

Elements in the configuration file are described concisely in comments, here they are described more thoroughly.

WIP

## Loci information

Information about loci, that are subjects for analysis by WarpSTR, must be described in the config file. An example is described `example/config.yaml`. Each loci must be defined by name and genomic coordinates. Then, you can either specify repeating motifs occuring in the locus in `motif` element, from which the input sequence for WarpSTR state automata is automatically created(this is recommended for starting users). The second way is to configure the input sequence by yourself in `sequence` element of the locus, however this is not a trivial task, so it is recommended for more advanced users. The other possibility is to use automatic configuration and then modify it by hand.
