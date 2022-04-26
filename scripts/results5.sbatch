#!/bin/bash

#SBATCH --account=NN8050K
#SBATCH --job-name=RES5
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=17
#SBATCH --mem-per-cpu=4G
#SBATCH --time=25:00:00

set -o errexit # Make bash exit on any error
set -o nounset # Treat unset variables as errors

## module restore
module --quiet purge
module load R/4.1.2-foss-2021b
## module load R/4.1.0-foss-2021a
## module load intel/2020b

cd ~/binclassfound

Rscript mcmc_5.R ${1} > mcmc_5_${1}.Rout 2>&1
