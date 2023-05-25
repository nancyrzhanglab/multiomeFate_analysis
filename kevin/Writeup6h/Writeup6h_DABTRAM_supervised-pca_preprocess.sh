#!/bin/bash
#$ -N spca_DABTRAM
#$ -j y
#$ -o ../../../../out/kevin/Writeup6h/qsub/
#$ -l m_mem_free=150G

Rscript --no-save Writeup6h_DABTRAM_supervised-pca_preprocess.R