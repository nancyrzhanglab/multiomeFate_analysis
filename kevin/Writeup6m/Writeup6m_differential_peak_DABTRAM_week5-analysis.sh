#!/bin/bash
#$ -N diff-peak_dabtram_w5-only
#$ -j y
#$ -o ../../../../out/kevin/Writeup6m/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup6m_differential_peak_DABTRAM_week5-analysis.R