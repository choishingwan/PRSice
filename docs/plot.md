# Introduction
PRSice will automatically generate the bar-plot and high-resolution plot
(if `--fastscore`) is not set. Quantile plots will also be generated if 
the `--quantile` parameter is provided

To see some example of PRSice output, please refer to [here](step_by_step.md#output-of-results)

!!! note
    These parameters, is not recognized by the PRSice binary

# Command

- `--bar-col-high [colour-code]`

   Colour of the most predicting threshold. Can either be a colour code like `red` or `\#E55738`, _e.g._ `--bar-col-high \#E55738`. Default: `firebrick`

- `--bar-col-lower [colour-code]`

  Colour of the poorest predicting threshold. Can either be a colour code like `blue` or `\#1290D9`, _e.g._ `--bar-col-high blue`. Default: `dodgerblue`

- `--bar-col-p`

  When set, will change the colour of bar to p-value threshold instead of
  the p-value from the association with phenotype

- `--bar-palatte`

  Colour palatte to be used for bar plotting when `--bar_col_p` is set. Default: `YlOrRd`

- `--multi-plot [N]`

  Plot the top N target phenotypes / gene sets in a summary plot, _e.g._ `--multi-plot 5`.


- `--plot`

  When set, will only perform plotting using existing PRSice result files. 
  All other parameters are still required such that PRSice
  can correctly locate the required input files for plotting.

- `--plot-set`

  The default behaviour of PRSet is to plot the bar-chart, high-resolution plot and
  quantile plot of the "Base" gene set, which consider
  all SNPs within the genome. By using the `--plot-set` option, you can plot the
  specific set of interest.


- `--quantile [N]` | `-q [N]`

    Number N of quantiles to plot, _e.g._ `--quantile 10`.
    No quantile plot will be generated when this is not provided.

-  `--quant-break`

    Parameter to indicate an uneven distribution of quantile. Values represent
    the upperbound of each quantile group. 

    e.g. With `--quantile 10 --quant-break 1,5,10`, the quantiles will be grouped into

    > $0\lt Q \le 1$, $1\lt Q \le 5$, $5\lt Q \le 10$


- `--quant-extract` | `-e`

  File containing sample ID to be plot on a separated
  quantile e.g. extra quantile containing only
  schizophrenia samples. Must contain IID. Should
  contain FID if `--ignore-fid` isn't set.

!!! note

    This will only work if the base and target has a different
    phenotype or if the target phenotype is quantitative


- `--quant-ref`

  Reference quantile for quantile plot. Default is number of quantiles divided by 2

  Or in the event where `--quant-break` is used, represent the upper bound of the 
  reference quantile

- `--scatter-r2`

  When set, will change the y-axis of the high resolution scatter plot to R2 instead
