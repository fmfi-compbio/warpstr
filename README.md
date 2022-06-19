# WarpSTR
WarpSTR is an alignment-free algorithm for analysing STR alleles using nanopore sequencing raw reads. The method uses guppy basecalling annotation output for the extraction of region of interest, and dynamic time warping based state automata for calling raw signal data. The tool can be configured to account for complex motifs containing interruptions and other variations such as substitutions or ambiguous bases. 

### Installation
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
```
python WarpSTR.py example/config.yaml
```

### Input data
Required input data are .fast5 files and .bam mapping files. In configuration file, the user is required to provide the path to the upper level path, in the `inputs` element. WarpSTR presumes that your data can come from multiple sequencing runs, but are of the same sample, and thus are to be analyzed together. For example, you have main directory for sample `subjectXY` with many subdirectories denoting sequencing runs i.e. `run_1`, `run_2`, with each run directory having its own .bam mapping file and .fast5 files. It is also possible to denote another path to input, in case of having data stored somewhere else (i.e. on the other mounted directory, as ONT data are very large), for example with the data from another run, i.e. `run_3`.

For the above example, `inputs` in the config could be defined as follows:
```
inputs:                       
  - path: /data/subjectXY 
    runs: run_1,run_2
  - path: /mnt/subjectXY 
    runs: run_3
```

Each directory as given by `path` and `runs`, i.e. `/data/subjectXY/run_1` and so on, is traversed by WarpSTR to find .bam files and .fast5 files.

## Output

... TBA

### Additional information
Newer .fast5 files are usually VBZ compressed, therefore VBZ plugin for HD5 is required to be installed, so WarpSTR can handle such files. See `https://github.com/nanoporetech/vbz_compression`. 