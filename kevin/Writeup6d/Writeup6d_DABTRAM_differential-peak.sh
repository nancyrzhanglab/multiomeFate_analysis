#!/bin/bash
#$ -N dabtram_peaks
#$ -j y
#$ -o ../../../../out/kevin/Writeup6d/qsub/
#$ -l m_mem_free=150G

Rscript --no-save Writeup6d_DABTRAM_differential-peak.R