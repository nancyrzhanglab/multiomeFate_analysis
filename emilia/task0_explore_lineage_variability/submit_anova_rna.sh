#!/bin/bash
#$ -N anova_rna
#$ -l m_mem_free=50G
#$ -j y

module load R
Rscript --vanilla anova_rna.R day10_COCL2_processed 