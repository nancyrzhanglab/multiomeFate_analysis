#!/bin/bash
#$ -N dabtram_d0_imputation-stepdown
#$ -j y
#$ -o ../../../../out/kevin/Writeup6p/qsub/
#$ -l m_mem_free=50G

Rscript --no-save Writeup6p_DABTRAM_day0_lineage-imputation_stepdown.R