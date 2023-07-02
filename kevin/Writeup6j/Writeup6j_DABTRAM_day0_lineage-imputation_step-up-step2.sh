#!/bin/bash
#$ -N DABTRAM_d0_stepup-step2
#$ -j y
#$ -o ../../../../out/kevin/Writeup6j/qsub/
#$ -l m_mem_free=25G

Rscript --no-save Writeup6j_DABTRAM_day0_lineage-imputation_step-up-step2.R