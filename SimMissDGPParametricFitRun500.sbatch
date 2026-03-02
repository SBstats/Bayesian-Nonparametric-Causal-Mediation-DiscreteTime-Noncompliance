#!/bin/bash
#SBATCH --job-name=EDPNoncomplianceProjParamSim                          # Job name
#SBATCH --ntasks=1
#SBATCH --mem=64GB
#SBATCH --mail-type=ALL                                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=XXXX                  # Where to send mail
#SBATCH --time=25:00:00                                 # Time limit hrs:min:sec
#SBATCH --array=1-500
#SBATCH --account=XXXX			
#SBATCH --qos=XXXX				
date; hostname; pwd

mkdir -p param_simu_results

# Load R module (module system is auto-sourced on compute nodes)
module load R

# Set library path for compute node packages (separate from interactive node)
export R_LIBS_USER=/home/XXXX/R/library_compute

Rscript Simulation-MissDGPParametricFitRun500.R