#!/bin/bash
#$ -N dabtram_d10_onlyrna
#$ -j y
#$ -o ../../../../out/kevin/Writeup6q/qsub/
#$ -l m_mem_free=20G

Rscript --no-save Writeup6q_DABTRAM-onlyRNA_day10_lineage-imputation.R