# Reproducibility for analyzing discretely-observed whale dives

This code repository contains data and workflows required to analyze whale dives
collected via satellite tags, which are often configured to collect
discretized depth measurements once every five minutes.

The repository contains simulation and analysis code, exploratory analyses, and
additional files used to manage analyses on a SLURM-based cluster.  The rest of
this document outlines how to reproduce the principal analyses.


# Developmental R packages required

The analysis uses two `R` packages that are not available on
[CRAN](https://cran.r-project.org), but can be installed with the help of the
[devtools](https://cran.r-project.org/web/packages/devtools/index.html) package.

```r
# Install package to assemble configuration files for analyses, plots, etc.
devtools::install_github("jmhewitt/composr")

# Install package to analyze dives on discrete spaces
devtools::install_github("jmhewitt/dsdive")
```


# All R packages required

You may need to install additional packages in order to run the scripts in this
repository.  The following code snippet uses the
[checkpoint](https://cran.r-project.org/web/packages/checkpoint/index.html)
package to scan the scripts in this repository to create a list of all
packages required.  The snippet then identifies and downloads packages that are
not already installed.

```r
library(checkpoint)

# scan for packages used in analysis
pkgs = scanForPackages()$pkg

# determine which packages are already installed
installed = pkgs %in% installed.packages()[,'Package']

# install missing packages
install.packages(pkgs[!installed])
```


# Steps for principal analyses

## Simulation study

  1. Run `scripts/sim_gen.R` to generate simulation data.
  2. Run `scripts/fit.R` to estimate model parameters for each simulation
     dataset.
     - `fit.R` is designed to use a configuration script to load a single
       dataset, specify prior distributions, and other sampler/output
       settings.  The configuration can be passed in as a command-line
       argument, or as the list variable `groups` in an interactive `R`
       session.  The simulations should use the configuration below, where
       `XXX` is replaced with the simulation dataset names
       `sim_tyack_more_known_end_30`, `sim_tyack_more_known_end_60`, and
       `sim_tyack_more_known_end_300`.

```r
groups = list(
  data = 'XXX',
  observation_model = 'exact_systematic',
  priors = 'tyack_simulation_priors',
  sampler = 'prod',
  subset = 'all_dives',
  validation= 'no_validation'
)
```


## Data analysis

```r
groups = list(
  data = 'zc84_800',
  observation_model = 'uniform_systematic',
  priors = 'tyack_priors_fixed_stage',
  sampler = 'test',
  subset = 'all_dives',
  validation= 'holdout_half'
)
```
