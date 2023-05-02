#!/bin/bash
#$ -N bino_extr_dabtram
#$ -j y
#$ -o ../../../../out/kevin/Writeup6g/qsub/
#$ -l m_mem_free=60G

Rscript --no-save Writeup6g_day0_binomial_extract-all_DABTRAM.R