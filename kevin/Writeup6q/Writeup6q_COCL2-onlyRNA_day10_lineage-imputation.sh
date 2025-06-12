#!/bin/bash
#$ -N cocl2_d10_onlyrna
#$ -j y
#$ -o ../../../../out/kevin/Writeup6q/qsub/
#$ -l m_mem_free=20G

Rscript --no-save Writeup6q_COCL2-onlyRNA_day10_lineage-imputation.R