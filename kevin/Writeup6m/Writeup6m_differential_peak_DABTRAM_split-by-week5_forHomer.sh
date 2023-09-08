#!/bin/bash
#$ -N diff-peak_w5dabtram
#$ -j y
#$ -o ../../../../out/kevin/Writeup6m/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup6m_differential_peak_DABTRAM_split-by-week5_forHomer.R