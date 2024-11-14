#!/bin/bash
#SBATCH --job-name=CHAMP_IMPORT
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem-per-cpu=64gb
#SBATCH --output=/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/logs/logs_CHAMP_IMPORT.out
#SBATCH --error=/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/logs/logs_CHAMP_IMPORT.err

module load gcc/12.1.0
module load miniconda3/23.1.0
source activate ChAMP
export TMPDIR=/gpfs/data/phs/groups/Projects/GEMS/james.li/conda_env/tmp

Rscript /gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/code/rscripts/CHAMP_IMPORT.R
