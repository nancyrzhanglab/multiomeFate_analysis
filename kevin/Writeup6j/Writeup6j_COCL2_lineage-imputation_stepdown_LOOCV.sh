#!/bin/bash
#$ -N COCL2_loocv
#$ -j y
#$ -o ../../../../out/kevin/Writeup6j/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup6j_COCL2_lineage-imputation_stepdown_LOOCV.R