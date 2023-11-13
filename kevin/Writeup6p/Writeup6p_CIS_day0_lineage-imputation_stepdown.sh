#!/bin/bash
#$ -N cis_d0_imputation-stepdown
#$ -j y
#$ -o ../../../../out/kevin/Writeup6p/qsub/
#$ -l m_mem_free=50G

Rscript --no-save Writeup6p_CIS_day0_lineage-imputation_stepdown.R