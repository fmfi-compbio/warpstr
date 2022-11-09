# Setting up input for WarpSTR

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
