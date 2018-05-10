# Introduction
PRSice will automatically generate the bar-plot and high-resolution plot
(if `--fastscore`) is not set. Quantile plots will also be generated if you
used the `--quantile` parameter.

To see some example of PRSice output, you can refer to [here](step_by_step.md#output-of-results)

!!! note
    These parameters, is not recognized by the PRSice binary

# Command

- `--bar-col-high`

   Colour of the most predicting threshold. Default: `firebrick`

- `--bar-col-lower`

  Colour of the poorest predicting threshold. Default: `dodgerblue`

- `--bar-col-p`

  When set, will change the colour of bar to p-value threshold instead of
  the p-value from the association with phenotype

- `--bar-palatte`

  Colour palatte to be used for bar plotting when `--bar_col_p` is set. Default: `YlOrRd`

- `--multi-plot`

  Plot the top N target phenotypes / gene sets in a summary plot


- `--plot`

  When set, will only perform plotting using existing PRSice result files.
  Users will still need to provide all other parameters such that PRSice
  can correctly locate the required input files for plotting.

- `--plot-set`

  The default behaviour of PRSet is to plot the bar-chart, high-resolution plot and
  quantile plot of the "Base" gene set, which consider
  all SNPs within the genome. By using the `--plot-set` option, you can plot the
  specific set of interest.


- `--quantile` | `-q`

    Number of quantiles to plot.
    No quantile plot will be generated when this is not provided.

- `--quant-extract` | `-e`

  File containing sample ID to be plot on a separated
  quantile e.g. extra quantile containing only
  schizophrenia samples. Must contain IID. Should
  contain FID if `--ignore-fid` isn't set.

!!! note

    This will only work if the base and target has a different
    phenotype f if the target phenotype is quantitative


- `--quant-ref`

  Reference quantile for quantile plot. Default is number of quantiles divided by 2

- `--scatter-r2`

  When set, will change the y-axis of the high resolution scatter plot to R2 instead
