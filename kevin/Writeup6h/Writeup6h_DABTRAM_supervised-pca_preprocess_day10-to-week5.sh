#!/bin/bash
#$ -N spca_DABTRAM_d10-w5
#$ -j y
#$ -o ../../../../out/kevin/Writeup6h/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup6h_DABTRAM_supervised-pca_preprocess_day10-to-week5.R