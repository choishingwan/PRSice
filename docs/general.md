# Introduction
PRSice is seperated into two main parts:

1. **PRSice executable**

    Responsible for the core algorithm


2. **PRSice Rscript**

    Responsible for generating the bar-plot, high-resolution plot and quantile plot

!!! tip

    You can perform PRS analysis using only the PRSice executable file
    (without the Rscript) with all non-plotting / R related parameters.
    This allow you to generate the PRS without performing the plotting,
    which might be useful for users who are interested in only obtaining
    the PRS.

## General Usage
```
usage: Rscript PRSice.R [options] <-b base_file> <-t target_file> <--prsice prsice_location>
```

## PRSice Binary (R specific)
- `--prsice`

    inform the location of the PRSice executable.
