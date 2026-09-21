#!/bin/bash
#SBATCH --mem=12G
#SBATCH - cwd
#SBATCH --output=mkfastq_out
#SBATCH --ntasks-per-node=16

module load cellranger/7.2.0-dntehee
module load bcl2fastq2/2.20.0.422-orocbiu

cd /data/lvd_qve/Projects/PTD_StrandSpecificCounting_scRNAseq/CellRanger
name=QVEU0056
path=$1
echo $path

cellranger mkfastq --run=${path} --csv=/data/lvd_qve/Projects/PTD_StrandSpecificCounting_scRNAseq/CellRanger/SampleSheet.csv --lanes=1,2 --localcores=16 --localmem=160 2>${name}_run.err 1>${name}_run.log
