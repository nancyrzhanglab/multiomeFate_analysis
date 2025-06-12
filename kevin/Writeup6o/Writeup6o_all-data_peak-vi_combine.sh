#!/bin/bash
#$ -N atac_peakvi_combine
#$ -j y
#$ -o ../../../../out/kevin/Writeup6o/qsub/
#$ -l m_mem_free=50G

Rscript --no-save Writeup6o_all-data_peak-vi_combine.R
