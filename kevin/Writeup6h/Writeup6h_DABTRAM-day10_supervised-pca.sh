#!/bin/bash
#$ -N spca_DABTRAM
#$ -j y
#$ -o ../../../../out/kevin/Writeup6h/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup6h_DABTRAM-day10_supervised-pca.R