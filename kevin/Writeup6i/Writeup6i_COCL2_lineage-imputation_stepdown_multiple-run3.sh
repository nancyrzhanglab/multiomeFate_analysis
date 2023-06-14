#!/bin/bash
#$ -N COCL2_imp_stpdwn3
#$ -j y
#$ -o ../../../../out/kevin/Writeup6i/qsub/
#$ -l m_mem_free=25G

Rscript --no-save Writeup6i_COCL2_lineage-imputation_stepdown_multiple-run3.R