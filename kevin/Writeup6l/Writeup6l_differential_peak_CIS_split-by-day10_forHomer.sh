#!/bin/bash
#$ -N diff-peak_cis_day10
#$ -j y
#$ -o ../../../../out/kevin/Writeup6l/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup6l_differential_peak_CIS_split-by-day10_forHomer.R