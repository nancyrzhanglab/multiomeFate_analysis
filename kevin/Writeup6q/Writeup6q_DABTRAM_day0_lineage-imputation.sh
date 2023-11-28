#!/bin/bash
#$ -N dabtram_d0_imputation
#$ -j y
#$ -o ../../../../out/kevin/Writeup6q/qsub/
#$ -l m_mem_free=50G

Rscript --no-save Writeup6q_DABTRAM_day0_lineage-imputation.R