# Introduction
PRSice is seperated into two main parts:

1. PRSice executable

    Responsible for the core algorithm


2. PRSice Rscript

    Responsible for generating the bar-plot, high-resolution plot and quantile plot

## General Usage
```
usage: Rscript PRSice.R [options] <-b base_file> <-t target_file> <--prsice prsice_location>
```

## PRSice Binary
- `--prsice`

    inform the location of the PRSice executable.

!!! tip
    You can also perform PRS using only the PRSice executable file
    (without the Rscript) with all non-plotting / R related parameters.
