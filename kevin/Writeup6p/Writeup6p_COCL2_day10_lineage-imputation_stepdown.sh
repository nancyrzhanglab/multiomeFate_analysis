#!/bin/bash
#$ -N cocl2_d10_imputation-stepdown
#$ -j y
#$ -o ../../../../out/kevin/Writeup6p/qsub/
#$ -l m_mem_free=50G

Rscript --no-save Writeup6p_COCL2_day10_lineage-imputation_stepdown.R