#!/bin/bash
#$ -N imputation_stepdown2
#$ -j y
#$ -o ../../../../out/kevin/Writeup6i/qsub/
#$ -l m_mem_free=20G

Rscript --no-save Writeup6i_DABTRAM_lineage-imputation_stepdown_multiple-run2.R