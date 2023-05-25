#!/bin/bash
#$ -N spca_DABTRAM_d0-d10
#$ -j y
#$ -o ../../../../out/kevin/Writeup6h/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup6h_DABTRAM_supervised-pca_preprocess_day0-to-day10.R