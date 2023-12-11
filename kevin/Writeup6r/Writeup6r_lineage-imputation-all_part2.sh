#!/bin/bash
#$ -N all_part2
#$ -j y
#$ -o ../../../../out/kevin/Writeup6r/qsub/
#$ -l m_mem_free=50G

Rscript --no-save Writeup6r_lineage-imputation-all_part2.R