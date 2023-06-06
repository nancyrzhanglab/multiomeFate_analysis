#!/bin/bash
#$ -N bin_extract_DABTRAM
#$ -j y
#$ -o ../../../../out/kevin/Writeup6i/qsub/
#$ -l m_mem_free=150G

Rscript --no-save Writeup6i_day0_binomial_extract-all_DABTRAM.R