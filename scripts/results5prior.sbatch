#!/bin/bash

#SBATCH --account=NN8050K
#SBATCH --job-name=RES5p
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=17
#SBATCH --mem-per-cpu=2G
#SBATCH --time=00:10:00

set -o errexit # Make bash exit on any error
set -o nounset # Treat unset variables as errors

## module restore
module --quiet purge
module load R/4.1.2-foss-2021b
## module load R/4.1.0-foss-2021a
## module load intel/2020b

cd ~/binclassfound

Rscript mcmc_5prior.R ${1} > mcmc_5prior_${1}.Rout 2>&1
