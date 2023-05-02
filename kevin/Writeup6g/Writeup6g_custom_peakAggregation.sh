#!/bin/bash
#$ -N peakAggregation
#$ -j y
#$ -o ../../../../out/kevin/Writeup6g/qsub/
#$ -l m_mem_free=150G

Rscript --no-save Writeup6g_custom_peakAggregation.R