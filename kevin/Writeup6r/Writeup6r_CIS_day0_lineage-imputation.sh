#!/bin/bash
#$ -N cis_d0
#$ -j y
#$ -o ../../../../out/kevin/Writeup6r/qsub/
#$ -l m_mem_free=20G

Rscript --no-save Writeup6r_CIS_day0_lineage-imputation.R