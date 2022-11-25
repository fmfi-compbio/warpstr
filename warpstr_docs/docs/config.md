# Configuration file

The only required parameter for running WarpSTR is the path to the config file. The file must be in YAML format. It must be populated with elements such as `inputs`, `output` and `reference_path`. An example is provided in `example/config.yaml`.

There are also many advanced parameters that are optional to set. List of all parameters are found in `example/advanced_params.yaml`. To set values for those parameters, just add those parameters to your main config and set them to the desired value. In other case, default values for those parameters are taken.

## Description of required parameters

Fields of the configuration file are described concisely in comments, here some of them are described more thoroughly.

### Input

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

!!!tip
    See [Test case](installation.md#test-case) for an example of test data and associated config.

### Reference path

Path to the reference fasta file. The reference is required for two things:

- in the `exp_signal_generation` step. Here the reference is used to extract flanking sequences of the locus.
- when `sequence` field is missing in locus configuration. Here, the reference is used to obtain the reference locus sequence and to automatically generate `sequence` field required for subsequent steps.

In other cases (skipping `exp_signal_generation` step and supplying `sequence` for all loci), the reference path is not needed.

### guppy_config

This field is a list, such as:

```yaml
guppy_config:
    path: guppy_executable
    flowcell: FLO-MIN106
    kit: SQK-LSK109
```

In `path` subfield supply the path to the executable guppy basecaller, and in `flowcell` and `kit` provide the nanopore sequencing information for the data.

This field must be filled in case of running step `guppy_annotation`, where input .fast5 files are basecalled (or basecalled again) to obtain the mapping between the basecalled sequence and the signal data.

## Loci information

Information about loci, that are subjects for analysis by WarpSTR, must be described in the config file. An example is described `example/config.yaml`:

Each locus has the following information:

- `name` - denotes the name of the output subdirectory for this locus.
- `coord` - genomic coordinates of the locus.
- `motif` - used to obtain the sequence for automaton automatically. Should be a pattern or a comma separated list of patterns such as `AGC,CGC`, that are presumed to repeat in the locus. This parameter is recommended for beginners or for initial analysis of the locus. For advanced usage see `sequence` field.
- `sequence` - used for constructing state automaton directly (flanks are automatically added accordingly). Supplying a better tailored sequence in this field produces usually more accurate results. See [Setting up automaton sequence](#automaton-sequence) for more information. Warning - `motif` field is ignored if using `sequence`.
- `noting` - currently not used, serves only as some kind of note.
- `flank_length` - It denotes how many nucleotides are taken from the reference as flanks. If not set, the default value of 110 is used.

An example for Huntington's disease locus:

```yaml
name: HD
coord: chr4:3,074,878-3,074,967
noting: AGC[19]AAC[1]AGC[1]CGC[1]CAC[1]CGC[7]
motif: AGC,CGC                      # use this
sequence: (AGC)AACAGCCGCCAC(CGC)    # or this
```

### Automaton sequence

The sequence for automaton is a concatenated string of three substrings: left flank, repeat sequence and right flank. All sequences are taken from the reference, but the repeat sequence can be set manually, mainly in cases when you presume that there are differences between the reference locus and your subject locus such as mutations in the repeating pattern or novel repeating patterns.

The repeat sequence has a strict syntax:

- IUPAC codes - denoting patterns
- parentheses `(` and `)` - everything between them is taken as the repeat structure that must occur at least once.
- curly brackets `{` and `}` - everything between them is taken as optional repeat structure that can occur zero times.
- all substrings that are not inside parentheses or brackets must occur exactly once.

For example, `(AGC)AACAGCCGCCAC(CGC)`, means that there could be 1-N repeats of AGC, followed by nonrepeating part AACAGCCGCCAC, and ending with 1-N repeats of (CGC).

It is important to not provide a repeat sequence that is too lenient, i.e. providing many possible states as that would make the aligning very slow and probably not accurate. This is due to the nature of the nanopore sequencing data, as many patterns are too similar to each other and due to noise or time dilation, alignment would not be accurate. For example, you could model `(AGC)AACAGCCGCCAC(CGC)` with a more lenient version such as `(MRC)`, where M denotes A or C, and R denotes A or G, but this sequence models a vast amount of other possibilities. Therefore, you should try to set the sequence as strict as possible.

!!! warning
    Do not use WarpSTR for homopolymeric sequences, i.e. do not set the sequence to contain parts such as `(C)`.
