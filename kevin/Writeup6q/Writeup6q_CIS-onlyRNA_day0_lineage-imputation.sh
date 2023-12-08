#!/bin/bash
#$ -N cis_d0_onlyrna
#$ -j y
#$ -o ../../../../out/kevin/Writeup6q/qsub/
#$ -l m_mem_free=20G

Rscript --no-save Writeup6q_CIS-onlyRNA_day0_lineage-imputation.R