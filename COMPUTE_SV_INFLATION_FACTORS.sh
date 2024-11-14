#!/bin/bash
#SBATCH --job-name=COMPUTE_SV_INFLATION_FACTORS
#SBATCH --output=/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/logs/logs_COMPUTE_SV_INFLATION_FACTORS_%A.out
#SBATCH --error=/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/logs/logs_COMPUTE_SV_INFLATION_FACTORS_%A.err
#SBATCH --time=36:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=32gb
#SBATCH --partition=tier2q

i=${ARGS1}

module load gcc/12.1.0
module load miniconda3/23.1.0
source activate r_env
export TMPDIR=/scratch/jll1/tmp

Rscript /gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/code/rscripts/COMPUTE_SV_INFLATION_FACTORS.R ${i} 
