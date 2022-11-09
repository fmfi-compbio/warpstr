# Advanced parameters

Many parameters can be set in addition to the parameters mentioned in the [Config section](config.md#configuration-file). The list of all parameters with their default values is in the repository `example/advanced_params.yaml`.
To change values for those parameters, simply add those parameters to your main config and set them to the desired value. Else, default values for those parameters are used.

Further below, we describe some of these parameters from the use case point of view.

## Using multiple threads

Add `threads` field and provide the desired value:

```yaml
threads: 32
```

!!!note
    Multiple threads are used for `tr_region_extraction` and `tr_region_calling` steps.

## Using nanopore data different from R9.4

WarpSTR by default uses the template pore model for the R9.4 chemistry. In other cases, such as as using R10 chemistry or other, you must provide the corresponding pore model (with expected signal values for all possible kmers) by path:

```yaml
pore_model_path: template_r10.model
```

!!!warning
    Different pore models are supported but they were not tested whether they produce satisfactory results.

## Changing hyperparameters for tandem repeat calling

WarpSTR offers the following parameters:
`spike_removal` - Possible values are `None` - no removal, `median3` and `median5` for median filtering window size of 3 and 5, respectively, and `Brute` for median filtering of window 5, but only for outliers, where outliers are determined as signal values not between 250 and 1000. Spike removal occur before normalization. Default: `Brute`
`min_values_per_state` - defines the minimum number of signal values required to map to one state. Default: 4. This is not recommended to change.
`states_in_segment` - defines the number of states that are contained in one window during the adaptive relaxation phase. Default 6. Lower values are less robust.
`min_state_similarity` - defines the minimum similarity between repeat states to raise a warning. Default: 0.75.

All these parameters are subfields of `tr_calling_config` field.

Further you can modify rescaling (polishing) hyperparameters:

`reps_as_one` - defines whether to merge all signal points of the repeating state together or not for the polishing phase. Default: False.
`threshold` - Max distance between state value and expected value to use the state for rescaling. Default: 0.5. Larger values mean that incorrectly aligned state (as the distance is larger) is used to rescale the data, which can results in less accurate rescaling.
`max_std` - Max standard deviation for signal points of the state to use that state for rescaling. Default: 0.5. Larger values mean that incorrectly aligned state (as the standard deviation is larger meaning it is too noisy) is used to rescale the data, which can results in less accurate rescaling.
`method` - Defines how to calculate state value from mapped signal points. Allowed `mean` and `median`.

For example, to change `spike_removal` to `median3` and rescaling hyperparameter threshold to 1, provide the following in your config:

```yaml
tr_calling_config:
    spike_removal: median3
    rescaling:
        threshold: 1
```

Further, you can turn off some visualizations (see [Output chapter](output.md)) by extending the `tr_calling_config`:

```yaml hl_lines="3-6"
tr_calling_config:
    spike_removal: median3
    visualize_alignment: False # signal alignments
    visualize_phase: False     # allele lengths split per phase
    visualize_strand: False    # allele lengths split per strand
    visualize_cost: False      # scatterplot of lengths vs state-wise cost
    rescaling:
        threshold: 1
```

## Changing hyperparameters for summarizing

For the summarizing step of WarpSTR, the following parameters can be modified:

- `min_weight` - Defines the minimum weight for heterozygosity. When one mixture model component has lower weight than this, the allele is declared as homozygous. Default: 0.2. Lower values mean that it is more probable to take a low frequency component as allele.
- `std_filter` - Defines the aggressiveness for filtering out predictions. Predictions are filtered if they are out of range: mean +- 2*std_filter. Default: 2. Higher values mean more aggressive filtering.
- `visualize` - Whether to output violinplot of clusters (see see [Output chapter](output.md)).
- `msa` - Whether to run [muscle](https://www.drive5.com/muscle/) to get multiple sequence alignment of predictions and basecalls. Default: False

To increase the aggressiveness of filtering, we can increase the `std_filter` hyperparameter as follows:

```yaml
genotyping_config:
    std_filter: 2.5
```

## Changing hyperparameters for alignment

For the `tr_extraction` step, we can modify the parameters for local alignment of flanks with the basecalled sequence. The optional parameters are:

- `accuracy_factor` - Defines whether the alignment score is good enough. A flank is deemed as not found when the ratio of its alignment score and the length is lower than this value. Default: 1.15. The value depends extremely on the scoring parameters.
- `identity_factor` - Defines whether the alignment score is good enough. A flank is deemed as not found when its alignment has lower identity than this, i.e. number of matches compared to the length. Default: 0.85. The value depends extremely on the scoring parameters.
- `match_score` - Defines how much is added to the alignment score when nucleotides match. Default: 2. Should be a positive value.
- `mismatch_score` - Defines how much is added to the alignment score when nucleotides do not match. Default: -3. Should be a negative value or smaller than `match_score`.
- `gap_open_score` - Defines how much is added to the alignment score when a gap is opened. Default: -3. Should be a negative value.
- `gap_extend_score` - Defines how much is added to the alignment score when a gap is extended. Default: -3. Should be a negative value smaller or equal to `gap_open_score`.

To change these parameters, we recommend to change all parameters at once, as they are strong dependencies. The example of change:

```yaml
alignment:
    accuracy_factor: 1.15
    identity_factor: 0.85
    match_score: 2
    mismatch_score: -3
    gap_open_score: -3
    gap_extend_score: -3
```

!!!warning
    These settings strongly influences whether the repeat signal would be localized correctly and whether the read would be used for subsequent repeat number characterization.
