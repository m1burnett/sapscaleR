# sapscaleR

`sapscaleR` provides models for up-scaling sap flux density measurements to
whole-tree volumetric sap flow rates using classic circular trunk theory and
the irregular-trunk theory proposed by Burnett et al. (2026) in *Agricultural
and Forest Meteorology* (<https://doi.org/10.1016/j.agrformet.2026.111357>).


## Installation

Once the GitHub repository is public, install with:

```r
install.packages("remotes")
remotes::install_github("m1burnett/sapscaleR")
```

## Quick Example

```r
library(sapscaleR)

curve <- make_beta_curve(mu = 0.348, K = 3.66)
trunk <- make_trunk(convexity = 0.8, diameter = 40)

SF_corr(
  curve,
  trunk,
  eval_mode = "z",
  dsw = 10,
  h_rel = 0.1
)
```

See `vignettes/sapscaleR-examples.Rmd` for a complete worked example. A knitted
HTML version is included under `docs/articles/sapscaleR-examples.html`.

## Data

Example data are bundled in `inst/extdata/` and can be accessed after
installation with:

```r
system.file("extdata", "pisonia-data.csv", package = "sapscaleR")
system.file("extdata", "avg.sfd.pigr.csv", package = "sapscaleR")
```

## Citation

If you use `sapscaleR` in published research, please cite Burnett et al. (2026).
The article DOI is <https://doi.org/10.1016/j.agrformet.2026.111357>.
