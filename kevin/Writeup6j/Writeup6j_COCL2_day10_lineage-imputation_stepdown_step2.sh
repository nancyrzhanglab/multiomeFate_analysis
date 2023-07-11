#!/bin/bash
#$ -N COCL2_day10-stepdown
#$ -j y
#$ -o ../../../../out/kevin/Writeup6j/qsub/
#$ -l m_mem_free=25G

Rscript --no-save Writeup6j_COCL2_day10_lineage-imputation_stepdown_step2.R