#!/bin/bash
#$ -N binomial_extract_key
#$ -j y
#$ -o ../../../../out/kevin/Writeup6g/qsub/
#$ -l m_mem_free=75G

Rscript --no-save Writeup6g_day0_binomial_extract-key.R