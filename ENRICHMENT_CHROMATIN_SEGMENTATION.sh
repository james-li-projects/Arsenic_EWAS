#!/bin/bash
#SBATCH --job-name=ENRICHMENT_CHROMATIN_SEGMENTATION
#SBATCH --output=/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/logs/logs_ENRICHMENT_CHROMATIN_SEGMENTATION_%A.out
#SBATCH --error=/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/logs/logs_ENRICHMENT_CHROMATIN_SEGMENTATION_%A.err
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32gb
#SBATCH --partition=tier2q

i=${ARGS1}

module load gcc/12.1.0
module load miniconda3/23.1.0
source activate r_env
export TMPDIR=/scratch/jll1/tmp

Rscript /gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/code/rscripts/ENRICHMENT_CHROMATIN_SEGMENTATION.R ${i} 
