#!/bin/bash
#SBATCH --ntasks-per-node = 16
#SBATCH --mem=10G
  
module load cellranger/7.2.0-dntehee
  
cd /data/lvd_qve/Projects/PTD_StrandSpecificCounting_scRNAseq/CellRanger
cellranger aggr --id scRNAseq_aggr_${2} --csv $1 --localcores 16 --localmem 100
