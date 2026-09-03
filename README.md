# Supplemental Materials for Liang & Barr, "Improving power without increasing sample size"

This repository contains the source materials (code and data) for the paper [*Improving power without increasing sample size: Model-based and design-based strategies for controlling time-dependent nuisance variation in experiments*](manuscript/Liang_Barr_Improving-Power.pdf) by Jinghui Liang and Dale J. Barr.

This repository also contains [supplementary materials](supplementary-materials.pdf) for the paper.

We have also created a Docker image [`dalejbarr/improving-exp-power`](https://hub.docker.com/repository/docker/dalejbarr/improving-exp-power/) that contains the software environment needed to exactly reproduce our results, including our manuscript.

## Getting started

### Recommended approach: Docker container

The recommended way to run the simulations is to use the Docker image, because it contains an identical software environment to that used for the simulations reported in the manuscript. The code below assumes you are just testing with only 2 monte carlo runs, and that you have a local directory *`/home/user/sims`* where you want the output to be stored. Substitute your own directory path and desired number of simulation runs.

``` shell
docker container run --rm --volume /home/user/sims:/output dalejbarr/improving-exp-power:latest do-all-simulations.sh 2
```

The output should look like below. After it finishes it saves the results to an R binary (`.rds`) file which you can read into a session using `readRDS()`.

```
using psrsim version '0.0.2.9000'
2 Monte Carlo runs with seed 1317376301 and np=25, nl=2, nr=6, ns=1
  (001/108) nl=02 nr=06 ns=01 eta=4 errs=expDecay bs=ps k=06... 0.4538174 secs 
  (002/108) nl=02 nr=06 ns=03 eta=4 errs=expDecay bs=ps k=06... 0.2768281 secs 
  (003/108) nl=02 nr=06 ns=06 eta=4 errs=expDecay bs=ps k=06... 0.3891978 secs
  [ . . . ]
  (106/108) nl=08 nr=24 ns=01 eta=4 errs=mixed bs=ps k=51... 28.4711 secs 
  (107/108) nl=08 nr=24 ns=12 eta=4 errs=mixed bs=ps k=51... 28.70183 secs 
  (108/108) nl=08 nr=24 ns=24 eta=4 errs=mixed bs=ps k=51... 29.6223 secs 
finished! processing results
results saved to '/output/D9ee42967a48b_1317376301_00002_25_2_06_2026-09-03 17:05:08.954532.rds'
total time elapsed: 10.4757826328278 mins
```

Most of the lines of output are there to show the progress of the simulations. A line such as:

```
(003/108) nl=02 nr=06 ns=06 eta=4 errs=expDecay bs=ps k=06... 0.3891978 secs
```

Means that parameter setting 3 of 108 was completed in about 0.39 seconds, and at this settings there were 2 levels of the independent variable (`nl=02`), with six repetitions at each level (`nr=06`) and six subblocks for PSR (`ns=06`) with exponential decay (learning) type error structure, with basis type P-splines (`bs=ps`) and 6 basis function (`k=06`).

**WARNING: It would take a very, very long time to run all the simulations (5000 Monte Carlo runs at each parameter setting) on a single computer.** For instance, it took around 10 minutes just to run 2 Monte Carlo runs at each of the 108 parameter settings on a standard laptop computer. The recommended approach if you want to run a large set would be to run smaller batches in parallel, e.g., on a computing cluster, and then combine the resulting `.rds` files.

If instead of running all the simulations you just want to poke around inside the container and possibly step through the `do-all-simulations.R` script, you can start the container with:

``` shell
docker container run --interactive --rm --tty dalejbarr/improving-exp-power:latest bash
```

and then `cd /var/scripts` to get to the directory with the R script. You can activate the emacs editor by typing `emacs do-all-simulations.R` at the command prompt and browse the code. The container has [Emacs Speaks Statistics (ESS)](https://ess.r-project.org/) so that you can fire up a console window within emacs and run commands line by line.

### Alternative approach: Run outside the container

If you are unable to use the Docker container, you can always try running the R script. First, make sure you have installed the `psrsim` package. This package requires that the [GNU Scientific Library](https://www.gnu.org/software/gsl/) is installed on your system, along with tools for compiling C++ code.

``` shell
devtools::install("psrsim")
```

You can run the R script from a command line like so:

``` shell
Rscript do-all-simulations.R nmc output-dir
```

where `nmc` is the number of Monte Carlo simulations and `output-dir` is the directory where the output file should be stored.

### Dealing with simulation results (`.rds`) files

The main script `do-all-simulations.R` creates an R binary file which you can read into an R session using the `readRDS()` function (or similarly, using `readr::read_rds` from the tidyverse). 

For instance:

``` shell
library("tidyverse")

dat <- read_rds("D9ee42967a48b_1317376301_00002_25_2_06_2026-09-03 17:05:08.954532.rds")
```

The object `dat` will have the format: 

```
# A tibble: 108 × 10
      np    nl    nr    ns errstr     etastep targ_k p_aov     p_lmm     p_gam 
   <int> <int> <int> <int> <chr>        <int>  <int> <list>    <list>    <list>
 1    25     2     6     1 expDecay         4      6 <dbl [2]> <dbl [2]> <dbl> 
 2    25     2     6     3 expDecay         4      6 <dbl [2]> <dbl [2]> <dbl> 
 3    25     2     6     6 expDecay         4      6 <dbl [2]> <dbl [2]> <dbl> 
 4    25     2     6     1 pinkNoise        4      6 <dbl [2]> <dbl [2]> <dbl> 
 5    25     2     6     3 pinkNoise        4      6 <dbl [2]> <dbl [2]> <dbl> 
 6    25     2     6     6 pinkNoise        4      6 <dbl [2]> <dbl [2]> <dbl> 
 7    25     2     6     1 randomWalk       4      6 <dbl [2]> <dbl [2]> <dbl> 
 8    25     2     6     3 randomWalk       4      6 <dbl [2]> <dbl [2]> <dbl> 
 9    25     2     6     6 randomWalk       4      6 <dbl [2]> <dbl [2]> <dbl> 
10    25     2     6     1 mixed            4      6 <dbl [2]> <dbl [2]> <dbl> 
# ℹ 98 more rows
# ℹ Use `print(n = ...)` to see more rows
```

Each row has results for a single parameter setting, where:

- `np` : number of participants (fixed at 25)
- `nl` : number of levels
- `nr` : number of repetitions per level
- `ns` : number of subblocks (for PSR; 1 = baseline)
- `errstr` : error structure
- `etastep` : effect size increment (eta); actual eta value varies depending on data structure (see manuscript)
- `targ_k` : starting number of basis functions
- `p_aov` : p-values for ANOVA (not reported in the manuscript); this is a 'list column' structure (nested vector)
- `p_lmm` : p-values for PSR (or baseline if `ns=1`)
- `p_gam` : p-values for the Wiggly Intercept model


## Compiling the manuscript

To compile the report from scratch, use:

``` shell
docker container run --rm --volume /home/user/sims:/output dalejbarr/improving-exp-power:latest make-manuscript
```

In the above command, for *`/home/user/sims`* substitute your own file path where you want to the PDF output to appear on your system.
	
To compile the supplementary materials, use:

``` shell
docker container run --rm --volume /home/user/sims:/output dalejbarr/improving-exp-power:latest make-supplementary
```

## Index of this repository

- `csl/`: APA7 citation style files needed for compiling the manuscript
- `Dockerfile`: the Docker file for building the container image. Available in case the container image ever disappears from Docker Hub.
- `docker_*_.sh`: Bash scripts that are copied into the container image (see Dockerfile)
- `manuscript/` : Files for generating the manuscript, which was written using [GNU Emacs](https://www.gnu.org/software/emacs/) [org-mode](https://orgmode.org).
- `psrsim/` : The psrsim package for R, used for data generation; can be installed locally using `devtools::install("psrsim")` in R. Please note that this will only work on a Linux-based operating system; if you are using any other operating system please use the Docker container.
- `results/` : simulation results; main file here is `results-power.rds`, which is an R binary file (load using `readRDS("results-power.rds")`).
- `scripts/` : contains the main R simulation script `do-all-simulations.R` which should be run in batch mode (non-interactively) with two command arguments:
  1. number of monte carlo runs;
  2. name of output directory where simulation results files will be stored 
- `supplementary/` : a report of supplementary simulations described in the manuscript

![CC-BY License](https://mirrors.creativecommons.org/presskit/buttons/88x31/png/by.png)
