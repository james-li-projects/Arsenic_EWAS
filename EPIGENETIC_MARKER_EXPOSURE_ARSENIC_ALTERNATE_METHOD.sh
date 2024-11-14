#!/bin/bash
#SBATCH --job-name=EPIGENETIC_MARKER_EXPOSURE_ARSENIC_ALTERNATE_METHOD
#SBATCH --output=/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/logs/EPIGENETIC_MARKER_EXPOSURE_ARSENIC_ALTERNATE_METHOD_%A.out
#SBATCH --error=/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/logs/EPIGENETIC_MARKER_EXPOSURE_ARSENIC_ALTERNATE_METHOD_%A.err
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=256gb
#SBATCH --partition=tier2q

i=${ARGS1}

module load gcc/12.1.0
module load miniconda3/23.1.0
source activate r_env
export TMPDIR=/scratch/jll1/tmp

Rscript /gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/code/rscripts/EPIGENETIC_MARKER_EXPOSURE_ARSENIC_ALTERNATE_METHOD.R ${i}
