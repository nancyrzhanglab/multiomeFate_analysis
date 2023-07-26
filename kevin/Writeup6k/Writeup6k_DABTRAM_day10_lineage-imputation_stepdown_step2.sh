#!/bin/bash
#$ -N DABTRAM_d10_stepdown_step2
#$ -j y
#$ -o ../../../../out/kevin/Writeup6k/qsub/
#$ -l m_mem_free=25G

Rscript --no-save Writeup6k_DABTRAM_day10_lineage-imputation_stepdown_step2.R