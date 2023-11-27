#!/bin/bash
#$ -N dabtram_d10_imputation
#$ -j y
#$ -o ../../../../out/kevin/Writeup6q/qsub/
#$ -l m_mem_free=50G

Rscript --no-save Writeup6q_DABTRAM_day10_lineage-imputation.R