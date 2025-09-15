# iComBat

This repository is dedicated to the code availability for the paper: [**iComBat: An Incremental Framework for Batch Effect Correction in DNA Methylation Array Data**](https://doi.org/10.1016/j.csbj.2025.09.014).

## Description

iComBat extends the standard ComBat method by allowing new batches to be corrected without re-processing the already corrected data.

## Requirements

- R ≥ 4.2
- CRAN / Bioconductor packages: `matrixStats`

## License

This repository is licensed under the [Artistic License 2.0](LICENSE).

## Usage

### ComBatInitialize

```r
ComBatInitialize(
  dat,
  batch,
  mod = NULL,
  mean.only = FALSE
)
```

**Arguments**

`dat`
Genomic measure matrix (dimensions probe x sample) - for example, expression matrix

`batch`
Batch covariate (only one batch allowed)

`mod`
Model matrix for outcome of interest and other covariates besides batch

`mean.only`
(Optional) FALSE If TRUE ComBat only corrects the mean of the batch effect (no scale adjustment)

**Value**

A list containing the following components:
- `n.batch`: Number of batches
- `n.genes.full`: Total number of genes including zero-variance genes
- `keep.rows`: Indices of genes kept for adjustment
- `zero.rows`: Indices of genes with zero variance within batches
- `grand.mean`: Grand mean across batches
- `var.pooled`: Pooled variance
- `n.covariates.original`: Original number of covariates
- `n.covariates`: Number of covariates after removing confounded ones
- `mod.coef`: Regression coefficients for covariates
- `mean.only`: Whether mean-only correction was used

### ComBatIncremental

```r
ComBatIncremental(
  dat.new,
  new.batch,
  fit,
  mod.new = NULL,
  par.prior = TRUE,
  prior.plots = FALSE,
  intercept.col = 1
)
```

**Arguments**

`dat.new`
Genomic measure matrix (dimensions probe x sample) for the new batch

`new.batch`
Batch identifier for the new data (for future extension)

`fit`
Output from ComBatInitialize function

`mod.new`
Model matrix for outcome of interest and other covariates besides batch for the new samples

`par.prior`
(Optional) TRUE indicates parametric adjustments will be used, FALSE indicates non-parametric adjustments will be used

`prior.plots`
(Optional) TRUE give prior plots with black as a kernel estimate of the empirical batch effect density and red as the parametric

`intercept.col`
(Optional) Column index of the intercept in mod.new (default: 1)

**Value**

data A probe x sample genomic measure matrix, adjusted for batch effects.

## Citation

Tomo, Y. & Nakaki, R. (2025). iComBat: An Incremental Framework for Batch Effect Correction in DNA Methylation Array Data. *Computational and Structural Biotechnology Journal*, In Press. 

© Yui Tomo, Ryo Nakaki (2025).
