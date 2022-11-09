# WarpSTR

WarpSTR is an alignment-free algorithm for analysing STR alleles using nanopore sequencing raw reads. The method uses guppy basecalling annotation output for the extraction of region of interest, and dynamic time warping based state automata for calling raw signal data. The tool can be configured to account for complex motifs containing interruptions and other variations such as substitutions or ambiguous bases.

See our preprint at: <https://www.biorxiv.org/content/10.1101/2022.11.05.515275v1>

See below for some quick steps how to install and run WarpSTR, or refer to more detailed [documentation](https://fmfi-compbio.github.io/warpstr/).

## Installation

WarpSTR can be easily installed using conda environment, frozen in `conda_req.yaml`. The conda environment can be created as follows:

```bash
conda env create -f conda_req.yaml
```

After installation, it is required to activate conda environment:

```bash
conda activate warpstr
```

WarpSTR was tested in Ubuntu 20.04 OS.

## Running WarpSTR

Required step to do before running WarpSTR is to prepare config file and add loci information.

### Config file

The input configuration file must be populated with elements such as `inputs`, `output` and `reference_path`. An example is provided in `example/config.yaml`.

There are also many advanced parameters that are optional to set. List of all parameters are found in `example/advanced_params.yaml`. To set values for those parameters, just add those parameters to your main config and set them to the desired value. In other case, default values for those parameters are taken.

### Loci information

Information about loci, that are subjects for analysis by WarpSTR, must be described in the config file. An example is described `example/config.yaml`. Each loci must be defined by name and genomic coordinates. Then, you can either specify repeating motifs occuring in the locus in `motif` element, from which the input sequence for WarpSTR state automata is automatically created(this is recommended for starting users). The second way is to configure the input sequence by yourself in `sequence` element of the locus, however this is not a trivial task, so it is recommended for more advanced users. The other possibility is to use automatic configuration and then modify it by hand.

### Running

After creating configuration file, running WarpSTR is simple as it requires only the path to the config file:

```bash
python WarpSTR.py example/config.yaml
```

### Input data

Required input data are .fast5 files and .bam mapping files. In configuration file, the user is required to provide the path to the upper level path, in the `inputs` element. WarpSTR presumes that your data can come from multiple sequencing runs, but are of the same sample, and thus are to be analyzed together. For example, you have main directory for sample `subjectXY` with many subdirectories denoting sequencing runs i.e. `run_1`, `run_2`, with each run directory having its own .bam mapping file and .fast5 files. It is also possible to denote another path to input, in case of having data stored somewhere else (i.e. on the other mounted directory, as ONT data are very large), for example with the data from another run, i.e. `run_3`.

For the above example, `inputs` in the config could be defined as follows:

```yaml
inputs:                       
  - path: /data/subjectXY 
    runs: run_1,run_2
  - path: /mnt/subjectXY 
    runs: run_3
```

Each directory as given by `path` and `runs`, i.e. `/data/subjectXY/run_1` and so on, is traversed by WarpSTR to find .bam files and .fast5 files.

## Output

The upper path for output is given in the .yaml configuration file as `output` element. Outputs are separated for each locus as subdirectories of this upper path, where names of subdirectories are the same as the locus name.

The output structure for one locus is as follows:

```bash
alignments/         # contains alignments of template flanks with reads
expected_signals/   # contains template flanks as sequences and expected signals
fast5/              # signals extracted as encompasssing the locus, stored as signle .fast5 files
predictions/        # contains visualizations of automaton alignments and basecalled sequences (see below)
summaries/          # contains visualizations produced in the last summarizing phase (see below)
overview.csv        # .csv file with read information and output
```

Some output files are optional and can be controlled by the .yaml config file.

### Predictions

In the `predictions` directory of each locus there would be a large variety of outputted files in other subdirectories.

In **basecalls** subdirectory are output files related to basecalling, such as `all.fasta` containing basecalled sequences of all reads encompassing the locus as given by SAM/BAM, `basecalls_all.fasta` containing only reads in which flanks were found. This file is further split per strand into `basecalls_reverse.fasta` and `basecalls_template.fasta`. In case of running muscle for MSA - multiple sequence alignment (controlled by advanced_params config), there would be `msa_all.fasta` file with MSA. In case of running summarizing, there would be `group1.fasta` and `group2.fasta` files where would be basecalled sequences split into groups as summarized by the last step of WarpSTR. In such case MSA output would be also created only for basecalled sequences of each group.

In `complex_repeat_units.csv` file there is counter for each repeat structure of the complex STR locus. Each row denote a read, and in columns are counts for repeat structures.

In **sequences** subdirectory there is analogous information as in **basecalls** subdirectory, but the information is not produced from the basecalled sequences but from sequences as given by WarpSTR.

In **DTW_alignments** subdirectory there are visualized alignments of STR signal with automaton (in both stages). Visualizations are truncated to first 2000 values.

### Summaries

In the `summaries` directory of each locus there is a myriad of optional visualizations:

- alleles.svg - Summarized predictions of repeat lengths in 1 or 2 groups and for WarpSTR and basecall.
- collapsed_predictions.svg - Complex repeat structure counts, only for WarpSTR.
- collapsed_predictions_strand.svg - As above, but further split by strand.
- complex_genotypes.svg - Summarized complex repeat structure counts in 1 or 2 groups.
- predictions_cost.svg - Scatterplot of state-wise cost and allele lengths.
- predictions_phase.svg - Violinplots of repeat lengths in the first and second phase.
- predictions_strand.svg - Violinplots of repeat lengths as split by strand.

## Additional information

Newer .fast5 files are usually VBZ compressed, therefore VBZ plugin for HD5 is required to be installed, so WarpSTR can handle such files. See `https://github.com/nanoporetech/vbz_compression`. 
