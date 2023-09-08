#!/bin/bash
#$ -N diff-peak_cocl2_d0-vs-w5
#$ -j y
#$ -o ../../../../out/kevin/Writeup6m/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup6m_differential_peak_COCL2_day0-vs-week5.R