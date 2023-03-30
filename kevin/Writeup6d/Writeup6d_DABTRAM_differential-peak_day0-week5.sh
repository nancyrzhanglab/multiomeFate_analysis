#!/bin/bash
#$ -N dt_peaks_d0-w5
#$ -j y
#$ -o ../../../../out/kevin/Writeup6d/qsub/
#$ -l m_mem_free=150G

Rscript --no-save Writeup6d_DABTRAM_differential-peak_day0-week5.R