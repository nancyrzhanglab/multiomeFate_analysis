#!/bin/bash
#$ -N cis_peaks_d10-w5
#$ -j y
#$ -o ../../../../out/kevin/Writeup6d/qsub/
#$ -l m_mem_free=150G

Rscript --no-save Writeup6d_CIS_differential-peak_day10-week5.R