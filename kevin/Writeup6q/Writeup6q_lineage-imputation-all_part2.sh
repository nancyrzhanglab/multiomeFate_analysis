#!/bin/bash
#$ -N all_imputation_part2
#$ -j y
#$ -o ../../../../out/kevin/Writeup6q/qsub/
#$ -l m_mem_free=25G

Rscript --no-save Writeup6q_lineage-imputation-all_part2.R