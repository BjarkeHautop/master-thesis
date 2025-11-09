# Master Thesis Scripts

This repository contains the code used to reproduce the results from my
Master's thesis.
Implementations of any example or analysis using Particle Marginal 
Metropolis-Hastings (PMMH) include precomputed results stored in `.rds` files. 
You can load these files directly without rerunning the simulations.

The PMMH code is built on my R package
[bayesSSM](https://bjarkehautop.github.io/bayesSSM/), which I developed 
alongside this thesis.  

A `renv` lock file is included to ensure consistent package versions. 
Additionally, a Dockerfile is provided to build a container and run the scripts
in a controlled environment.

## Running the R scripts locally

1. Clone the repository:

```powershell
git clone https://github.com/BjarkeHautop/master-thesis
cd master-thesis/R
```

2. Restore the `renv` environment:

```r
# install.packages("renv")
renv::restore()
```

3. Run a specific script (e.g., `Example 3.1/Example 3.1.R`) in `R`:

```r
source("Example 3.1/Example 3.1.R")
```

Replace `Example 3.1/Example 3.1.R` with the path to any other script you 
wish to run.

## Running the R scripts using Docker

1. Clone the repository:

```powershell
git clone https://github.com/BjarkeHautop/master-thesis
cd master-thesis/R
```

2. Build the Docker image (if not already built):

```powershell
docker build -t master-thesis .
```

3. Run a specific script (e.g. `Example 3.1/Example 3.1.R`) like this:

```powershell
docker run --rm -v /full/path/to/repo/R:/project master-thesis Rscript "/project/Example 3.1/Example 3.1.R"
```

Replace the path with the actual path to the cloned repository on your
machine. Replace `Example 3.1/Example 3.1.R` with any other script you wish to 
run.
