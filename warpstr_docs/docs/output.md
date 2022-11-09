# Output

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

## Predictions

In the `predictions` directory of each locus there would be a large variety of outputted files in other subdirectories.

In **basecalls** subdirectory are output files related to basecalling, such as `all.fasta` containing basecalled sequences of all reads encompassing the locus as given by SAM/BAM, `basecalls_all.fasta` containing only reads in which flanks were found. This file is further split per strand into `basecalls_reverse.fasta` and `basecalls_template.fasta`. In case of running muscle for MSA - multiple sequence alignment (controlled by advanced_params config), there would be `msa_all.fasta` file with MSA. In case of running summarizing, there would be `group1.fasta` and `group2.fasta` files where would be basecalled sequences split into groups as summarized by the last step of WarpSTR. In such case MSA output would be also created only for basecalled sequences of each group.

In `complex_repeat_units.csv` file there is counter for each repeat structure of the complex STR locus. Each row denote a read, and in columns are counts for repeat structures.

In **sequences** subdirectory there is analogous information as in **basecalls** subdirectory, but the information is not produced from the basecalled sequences but from sequences as given by WarpSTR.

In **DTW_alignments** subdirectory there are visualized alignments of STR signal with automaton (in both stages). Visualizations are truncated to first 2000 values.

## Summaries

In the `summaries` directory of each locus there is a myriad of optional visualizations:

- alleles.svg - Summarized predictions of repeat lengths in 1 or 2 groups and for WarpSTR and basecall.
- collapsed_predictions.svg - Complex repeat structure counts, only for WarpSTR.
- collapsed_predictions_strand.svg - As above, but further split by strand.
- complex_genotypes.svg - Summarized complex repeat structure counts in 1 or 2 groups.
- predictions_cost.svg - Scatterplot of state-wise cost and allele lengths.
- predictions_phase.svg - Violinplots of repeat lengths in the first and second phase.
- predictions_strand.svg - Violinplots of repeat lengths as split by strand.