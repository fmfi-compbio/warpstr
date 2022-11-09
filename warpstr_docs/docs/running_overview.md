# Running WarpSTR

To run WarpSTR, first you must create a configuration file, that must adhere to YAML syntax. Check your input .fast5 and .bam files and populate the config accordingly (See [Config section](config.md)). Also, you can change some advanced parameters, for that see [Advanced parameters section](advanced_params.md)). Then, running WarpSTR is simple as it requires only the path to the config file. You can run WarpSTR for as many loci as you want, but all the information must be in the config.

See `example/config.yaml` for an example of config file. After populating the config file with your values, run WarpSTR as follows:

```bash
python WarpSTR.py example/config.yaml
```
