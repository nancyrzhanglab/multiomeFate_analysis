#!/bin/bash
#$ -N diff-peak_w5cocl2
#$ -j y
#$ -o ../../../../out/kevin/Writeup6m/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup6l_differential_peak_COCL2_split-by-day10_forHomer_part2.R