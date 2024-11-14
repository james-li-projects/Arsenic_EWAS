#!/bin/bash
#SBATCH --job-name=RUN_EWAS_HEALS_COMBINED_celltype_noSV_log2_UrAsgmCr_log2_WArsenic
#SBATCH --output=/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/logs/logs_RUN_EWAS_HEALS_COMBINED_celltype_noSV_log2_UrAsgmCr_log2_WArsenic_%A_%a.out
#SBATCH --error=/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/logs/logs_RUN_EWAS_HEALS_COMBINED_celltype_noSV_log2_UrAsgmCr_log2_WArsenic_%A_%a.err
#SBATCH --time=90:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=32gb

i=${ARGS1}

module load gcc/12.1.0
module load miniconda3/23.1.0
source activate test_empty
export TMPDIR=/scratch/jll1/tmp

Rscript /gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/code/rscripts/QSUB_EWAS_HEALS_COMBINED_celltype_noSV_log2_UrAsgmCr_log2_WArsenic.R ${i} 