#!/bin/bash
#$ -N cocl2_d0_imputation-stepdown_step2
#$ -j y
#$ -o ../../../../out/kevin/Writeup6n/qsub/
#$ -l m_mem_free=50G

Rscript --no-save Writeup6n_COCL2_day0_lineage-imputation_stepdown_step2.R