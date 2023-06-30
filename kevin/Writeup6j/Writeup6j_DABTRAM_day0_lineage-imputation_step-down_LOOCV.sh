#!/bin/bash
#$ -N DABTRAM_d0_stepdown_loocv
#$ -j y
#$ -o ../../../../out/kevin/Writeup6j/qsub/
#$ -l m_mem_free=25G

Rscript --no-save Writeup6j_DABTRAM_day0_lineage-imputation_step-down_LOOCV.R