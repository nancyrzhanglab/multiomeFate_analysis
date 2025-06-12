#!/bin/bash
#$ -N corr_export
#$ -j y
#$ -o ../../../../out/kevin/Writeup6n/qsub/
#$ -l m_mem_free=25G

Rscript --no-save Writeup6n_lineage-imputation_day0-day10-export.R