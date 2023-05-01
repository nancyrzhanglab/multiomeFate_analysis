#!/bin/bash
#$ -N binomial_extract
#$ -j y
#$ -o ../../../../out/kevin/Writeup6g/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup6g_day0_binomial_extract-all.R