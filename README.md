# WarpSTR

WarpSTR is an alignment-free algorithm for analysing STR alleles using nanopore sequencing raw reads. The method uses guppy basecalling annotation output for the extraction of region of interest, and dynamic time warping based state automata for calling raw signal data. The tool can be configured to account for complex motifs containing interruptions and other variations such as substitutions or ambiguous bases.

See our preprint at: <https://www.biorxiv.org/content/10.1101/2022.11.05.515275v1>

See below for some quick steps how to install and run WarpSTR, or refer to more detailed [documentation](https://fmfi-compbio.github.io/warpstr/).

## 1 Installation

WarpSTR can be installed using conda or pipenv. To install conda, please follow [the official guide](https://conda.io/projects/conda/en/latest/user-guide/install/index.html). To install pipenv, simple `pip install pipenv` should suffice.

WarpSTR was tested in Ubuntu 20.04 and Ubuntu 22.04. Used Python version is 3.7.

### 1.a) Installing using conda

Clone this repository. Then, create the conda environment:

```bash
conda env create -f conda_req.yaml
```

After installation, it is required to activate conda environment:

```bash
conda activate warpstr
```

### 1.b) Installing using pipenv

Clone this repository. The pipenv environment can be installed from Pipfile.lock as follows:

```bash
pipenv sync
```

After installation, it is required to activate the environment:

```bash
pipenv shell
```

## 2 Running the test case

In `test/test_input` there is a small test dataset with 10 reads for one locus. There is also the template for config file required by WarpSTR, `test/config_template.yaml`. You can check whether WarpSTR works correctly simply by running:

```bash
bash run_test_case.sh
```

When running this wrapper script, the script will prompt you to provide the required paths and run the WarpSTR for you using the test data. Output files will be then stored in `test/test_output/` as given in the config file. The script should take approximately 3-5 minutes and at the end, you should see something like:

```text
Results stored in overview file XY
Allele lengths as given by WarpSTR: (44, 40)
```

## 3 Running WarpSTR

Running WarpSTR is simple as it requires only the path to the configuration file:

```bash
python WarpSTR.py example/config.yaml
```

WarpSTR consists of multiple complex steps doing the following:

1. extracting reads encompassing the locus coordinates - requires BAM mapping files and multi .fast5.
2. extracting STR regions from reads - requires Guppy basecaller.
3. determining the alelle length for reads.
4. genotyping alelle lengths and determining zygosity.

If you want to run a whole WarpSTR pipeline then continue reading, else skip to the [WarpSTR steps](#5-warpstr-steps).

### 3.1 Config file

In the input configuration file (see `example/config.yaml` for an example) you must set the following elements:

- `reference_path` - path to the fasta file - the reference genome, the same which was used for mapping basecalled reads.
- `guppy_config` - path to the executable Guppy basecaller and info about the sequencing (flowcell and kit).
- `output` - path to the directory, when WarpSTR will produce output results.
- `loci` - loci information, see [below](#32-loci-information).
- `inputs` - input data, see [below](#33-input-data).

There are also many advanced parameters that are optional to set. List of all parameters are found in `example/advanced_params.yaml`. To set values for those parameters, just copy the elements to your main config and change valeus to your desired values. In other case, default values for those parameters are taken.

### 3.2 Loci information

Information about loci, that are subjects for analysis by WarpSTR, must be described in the config file. An example is described `example/config.yaml`. Each locus must be defined by name and genomic coordinates (these must match with the reference), and either motif or sequence:

```yaml
name: HD
coord: chr4:3,074,878-3,074,967
motif: AGC,CGC
# sequence: (AGC)AACAGCCGCCAC(CGC)
```

The `motif` element is recommended for beginners, as the input sequence for WarpSTR state automata is automatically created. In this element, possible repeat units should be provided.

The second way is to configure the input sequence for automata by yourself in the `sequence` element of the locus. This is not a trivial task, so it is recommended for more advanced users. The other possibility is to use automatic configuration and then modify it by hand. See the preprint for the information about the state automata.

### 3.3 Input data

Required input data are .fast5 files and .bam mapping files. In configuration file, the user is required to provide the path to the upper level path, in the `inputs` element. WarpSTR presumes that your data can come from multiple sequencing runs, but are of the same sample, and thus should be analyzed together, see [documentation](https://fmfi-compbio.github.io/warpstr/) in that case.

The simple case is like in the test case:

```bash
test_input/
└── test_run1
    ├── fast5s
    │   └── batch_0.fast5
    └── mapping
        ├── mapping.bam
        └── mapping.bam.bai
```

The names `test_run1` and `test_input` are then used in the configuration file for the `inputs` element:

```yaml
inputs:                       
  - path: test/test_input
    runs: test_run1
```

Names of subdirectories such as `fast5s` and `mapping` are not important, but .fast5 files and .bam files must have the correct extension.

## 4 Output

The upper path for output is given in the .yaml configuration file as `output` element. Each locus has the separate output - a new subdirectory of this upper path with locus name is created, where the output is stored.

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

### 4.1 Predictions

In the `predictions` directory of each locus there would be a large variety of outputted files in other subdirectories.

In **basecalls** subdirectory are output files related to basecalling, such as `all.fasta` containing basecalled sequences of all reads encompassing the locus as given by SAM/BAM, `basecalls_all.fasta` containing only reads in which flanks were found. This file is further split per strand into `basecalls_reverse.fasta` and `basecalls_template.fasta`. In case of running muscle for MSA - multiple sequence alignment (controlled by advanced_params config), there would be `msa_all.fasta` file with MSA. In case of running summarizing, there would be `group1.fasta` and `group2.fasta` files where would be basecalled sequences split into groups as summarized by the last step of WarpSTR. In such case MSA output would be also created only for basecalled sequences of each group.

In `complex_repeat_units.csv` file there is counter for each repeat structure of the complex STR locus. Each row denote a read, and in columns are counts for repeat structures.

In **sequences** subdirectory there is analogous information as in **basecalls** subdirectory, but the information is not produced from the basecalled sequences but from sequences as given by WarpSTR.

In **DTW_alignments** subdirectory there are visualized alignments of STR signal with automaton (in both stages). Visualizations are truncated to first 2000 values.

### 4.2 Summaries

In the `summaries` directory of each locus there is a myriad of optional visualizations:

- alleles.svg - Summarized predictions of repeat lengths in 1 or 2 groups and for WarpSTR and basecall.
- collapsed_predictions.svg - Complex repeat structure counts, only for WarpSTR.
- collapsed_predictions_strand.svg - As above, but further split by strand.
- complex_genotypes.svg - Summarized complex repeat structure counts in 1 or 2 groups.
- predictions_cost.svg - Scatterplot of state-wise cost and allele lengths.
- predictions_phase.svg - Violinplots of repeat lengths in the first and second phase.
- predictions_strand.svg - Violinplots of repeat lengths as split by strand.

## 5 WarpSTR steps

WarpSTR pipeline steps are toggleable in the config file, i.e. you can skip them, by turning them to False:

```yaml
single_read_extraction: True   # Extracts reads mapped to the locus and stores them in single .fast5 format
guppy_annotation:       True   # Annotates .fast5 files with mapping between basecalled sequence and the signal
exp_signal_generation:  True   # Generates expected signals for flanks and repeats
tr_region_extraction:   True   # Finds tandem repeat region in read using alignment of basecalled sequence and reference repeat sequence
tr_region_calling:      True   # Uses state automata with DTW alignment to find the number of repeats for each signal
genotyping:             True   # Predicts the final allele lengths from the predicted repeat numbers of each read 
```

### 5.1 Extraction of locus reads

Here, .BAM and multi-fast5 files are required. The following config elements must be set:

- `inputs` element - defining directories containing .BAM and .fast5 files
- `loci` element - defining genomic coordinates
- `single_read_extraction` element set to `True`

In the output directory (given by `output` element) the state of the locus output subdirectory after running this step would be:

```tree
{locus_name}
├── fast5
│   └── {run_id}
│       ├── {read_name1}.fast5
│       ├── {read_name2}.fast5
│       └── ...
└── overview.csv - index of extracted reads with `name`,`run_id`,`reverse` values for each read
```

#### Skipping this step

If you have already single .fast5s ready for the locus and want to skip this step, you should simulate the outcome of the first step:

1. Create the subdirectory in the output directory with the same as the name of the locus in the config.
2. In the locus subdir create the `fast5/run_id` directory, where you copy single .fast5 reads (See above the output example)
3. In the locus subdir create `overview.csv` file where for each read signal there should be a row with three columns: `name`,`run_id`,`reverse`, Where `name` is the name as the read_name, and `reverse` having either True or False value, denoting the strand.

For example, the overview.csv for the above case would be:

```csv
read_name,run_id,reverse
read_name1,run_id,False
read_name2,run_id,True
...
```

Then, do not forget to turn off the step in the config file:

```yaml
single_read_extraction: False   # Extracts reads mapped to the locus and stores them in single .fast5 format
```

### 5.2 Extraction of STR regions

Requires executable Guppy basecaller (and completed previous pipeline step).

In this step, reads are basecalled again so they would be annotated with the mapping between basecalls and signal values. This mapping is then used to localize the STR region in signals.

The state of the locus output directory after running this step would be:

```tree
{locus_name}
├── fast5
│   └── {run_id}
│       ├── annot
│       │   ├── {read_name1}.fast5
│       │   ├── {read_name2}.fast5
│       │   └── ...
└── overview.csv - index of extracted reads with `name`,`run_id`,`reverse` values for each read.
In addition, there would be 'l_start_raw', 'r_end_raw' values, corresponding to signal positions, where starts the left flank and ends the right flank.
```

#### Skipping this step

If you have already .fast5 signals with localized STR regions, you again must simulate the output of this step. The other option is to use our script `prepare_caller_only.py`. It requires two things:

`config file` - the same as you would use further. The important thing is to set the `output` and `loci`
`.csv file` - with one row for .fast5 signal, and these required columns:

- `fast5_path` - path to the .fast5 read path.
- `locus` - name of the locus associated with the read.
- `read_name` - name of the read.
- `reverse` - strand information, True or False.
- `l_start_raw` - signal position where the left flank starts.
- `r_end_raw` - signal position where the right flank ends.

The example case is in the repository in `test/test_caller_only`. To run, provide the reference_path there in the config, and run using:

```bash
python prepare_caller_only.py --config test/test_caller_only/config_caller_only.yaml --file test/test_caller_only/example.csv
```

This creates a simulated output of the previous step in the `test/test_caller_only/Human_STR_1108232`. Then you can run WarpSTR:

```bash
python WarpSTR.py test/test_caller_only/config_caller_only.yaml
```

#### Important notes

- The `l_start_raw` and `r_end_raw` can be set approximately, 100-200 positions off should pose no problem for the correct result.
- The `l_start_raw` and `r_end_raw` must correspond to the flank positions, i.e. the flank length must be set to the same value in the config for `loci`.
- We currently do not support direct input of signal values of STR for STR calling.

## 6 Additional information

Newer .fast5 files are usually VBZ compressed, therefore VBZ plugin for HD5 is required to be installed, so WarpSTR can handle such files. See `https://github.com/nanoporetech/vbz_compression`.
