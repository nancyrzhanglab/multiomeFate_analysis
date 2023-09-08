#!/bin/bash
#$ -N diff-peak_cocl2_w5-only
#$ -j y
#$ -o ../../../../out/kevin/Writeup6m/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup6m_differential_peak_COCL2_week5-analysis.R