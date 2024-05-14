#!/bin/bash#
#SBATCH --ntasks-per-core=8
#SBATCH --mem=125G

module load cellranger/7.2.0-dntehee

cd /data/lvd_qve/Projects/PTD_StrandSpecificCounting_scRNAseq/CellRanger/
template=/data/lvd_qve/QVEU_Code/sequencing/template_fastas/refdata-gex-GRCh38-2020-A/fasta/GRCh38-2020-A_PV_GTF_mRuby/
runID=${1}

cellranger count --sample=Mock_5h --transcriptome=$template --id=Mock_5h --fastqs=/data/lvd_qve/Projects/PTD_StrandSpecificCounting_scRNAseq/81123_SCISSORS_FreshAnalysis/PV_Rep/CellRanger_Analysis/AAC2L5HM5_12/outs/fastq_path/ --localcores=16 --localmem=125 --output-dir /data/lvd_qve/Projects/PTD_StrandSpecificCounting_scRNAseq/CellRanger/Mock_5h_PV
