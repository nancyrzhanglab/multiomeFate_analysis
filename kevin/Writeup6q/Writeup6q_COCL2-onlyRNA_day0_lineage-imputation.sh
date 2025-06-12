#!/bin/bash
#$ -N cocl2_d0_onlyrna
#$ -j y
#$ -o ../../../../out/kevin/Writeup6q/qsub/
#$ -l m_mem_free=20G

Rscript --no-save Writeup6q_COCL2-onlyRNA_day0_lineage-imputation.R