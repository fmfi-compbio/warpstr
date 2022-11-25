# Installation

First, clone the repository:

```bash
git clone git@github.com:fmfi-compbio/warpstr.git
```

WarpSTR can then be easily installed using the conda environment, frozen in `conda_req.yaml`. The conda environment can be created as follows:

```bash
conda env create -f conda_req.yaml
```

After installation, it is required to activate conda environment:

```bash
conda activate warpstr
```

WarpSTR was tested in Ubuntu 20.04 OS.

If you have any problems with installation, please raise an issue on Github and we will look further into it as soon as possible.

## Test case

In `test/test_input` there is a small test dataset. There are 10 reads for one locus. In the template config `test/config_template.yaml` you can see that inputs are defined as follows:

```yaml
inputs:                       
  - path: test/test_input
    runs: test_run1
```

In the test data, there is only one run, called `test_run1` with all sequencing files stored there (i.e. BAM files and .fast5 files):

```bash
test_input/
└── test_run1
    ├── fast5s
    │   └── batch_0.fast5
    └── mapping
        ├── mapping.bam
        └── mapping.bam.bai
```

!!!note
    In WarpSTR, we allow for multiple sequencing runs, which are defined by `runs` element, where other runs are defined by comma, i.e. `test_run1,test_run2` and so on. WarpSTR concatenates `path` value with each comma splitted value of `runs` to obtain a list of paths, that are then search recursively for input files (.BAM and .fast5 files).

To run test data, you need to define `reference_path` and `guppy_config.path` in the config, or you can use provided wrapper bash script `run_test_case.sh`. Simply, when in WarpSTR directory (and with activated conda environment), run the wrapper script:

```bash
bash run_test_case.sh
```

When running this wrapper script, the script will prompt you to provide the required paths and run the WarpSTR for you. Output files will be then stored in `test/test_output/` as given in the config file.

!!!warning
    Test data were produced using human genome reference denoted as `GRCh38`, so ensure that you provide the path to the refernce of this version. You can download the reference for example from <https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.26/>
