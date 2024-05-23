#!/bin/sh

srun --pty bash
module load cellranger/7.2.0-dntehee
cd /data/lvd_qve/QVEU_Code/sequencing/template_fastas/refdata-gex-GRCh38-2020-A_PV_GFP_mRuby/
cellranger mkref --genome=GRCh38-2020-A_PV_GTF_mRuby --fasta=human_pv_mrubygfp.fa --genes=genomePVmrubygtf.gtf
