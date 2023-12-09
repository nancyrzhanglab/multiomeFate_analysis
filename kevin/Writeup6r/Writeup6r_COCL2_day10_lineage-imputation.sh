#!/bin/bash
#$ -N cocl2_d10
#$ -j y
#$ -o ../../../../out/kevin/Writeup6r/qsub/
#$ -l m_mem_free=20G

Rscript --no-save Writeup6r_COCL2_day10_lineage-imputation.R