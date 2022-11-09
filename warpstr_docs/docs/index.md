# WarpSTR

WarpSTR is an alignment-free algorithm for analysing STR alleles using nanopore sequencing raw reads. The method uses guppy basecalling annotation output for the extraction of region of interest, and dynamic time warping based state automata for calling raw signal data. The tool can be configured to account for complex motifs containing interruptions and other variations such as substitutions or ambiguous bases.

See our preprint at: <https://www.biorxiv.org/content/10.1101/2022.11.05.515275v1>

## Documentation overview

See [Installation chapter](installation.md) for information about how to install WarpSTR.

To know more about how to setup the configuration file, see [Setting up the config](config.md).

Before running, ensure that input files are located accordingly as set up in the config file. See [Input section](input.md).

To know more about output files and what they represent, see [Output chapter](output.md).
