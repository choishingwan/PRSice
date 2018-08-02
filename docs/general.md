# Introduction
PRSice is seperated into two main parts:

1. **PRSice executable**

    Responsible for the core algorithm


2. **PRSice Rscript**

    Responsible for generating the bar-plot, high-resolution plot and quantile plot

!!! tip

    PRS analysis can be performed using only the PRSice executable file
    (without the Rscript) with all non-plotting / R related parameters.
    This generates the PRS without plotting, which might be useful for 
    those who are only interested in obtaining the PRS.

## General Usage
```
usage: Rscript PRSice.R [options] \
    --base <base_file> \
    --target < target_file> \
    --prsice <prsice_location>
```

## PRSice Binary (R specific)
- `--prsice`

    Location of the PRSice executable.
