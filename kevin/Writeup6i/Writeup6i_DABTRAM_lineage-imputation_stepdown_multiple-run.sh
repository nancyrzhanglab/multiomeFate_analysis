#!/bin/bash
#$ -N imputation_stepdown
#$ -j y
#$ -o ../../../../out/kevin/Writeup6i/qsub/
#$ -l m_mem_free=50G

Rscript --no-save Writeup6i_DABTRAM_lineage-imputation_stepdown_multiple-run.R