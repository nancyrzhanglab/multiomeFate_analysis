#!/bin/bash
#$ -N COCL2_day10-stepup
#$ -j y
#$ -o ../../../../out/kevin/Writeup6j/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup6j_COCL2_day10_lineage-imputation_stepwise-up.R