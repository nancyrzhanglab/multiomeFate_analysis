#!/bin/bash
#$ -N DABTRAM_day0_extract-growth
#$ -j y
#$ -o ../../../../out/kevin/Writeup6k/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup6k_DABTRAM_day0_lineage-imputation_stepdown_LOOCV_postprocess-extract.R