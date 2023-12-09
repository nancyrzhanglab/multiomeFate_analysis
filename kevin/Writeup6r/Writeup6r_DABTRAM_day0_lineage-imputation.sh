#!/bin/bash
#$ -N dabtram_d0
#$ -j y
#$ -o ../../../../out/kevin/Writeup6r/qsub/
#$ -l m_mem_free=20G

Rscript --no-save Writeup6r_DABTRAM_day0_lineage-imputation.R