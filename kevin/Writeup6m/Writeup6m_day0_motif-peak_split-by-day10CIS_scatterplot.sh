#!/bin/bash
#$ -N day0_motif-peak_d10cis_scatterplot
#$ -j y
#$ -o ../../../../out/kevin/Writeup6m/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup6m_day0_motif-peak_split-by-day10CIS_scatterplot.R