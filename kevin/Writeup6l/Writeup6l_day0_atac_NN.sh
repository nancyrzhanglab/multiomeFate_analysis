#!/bin/bash
#$ -N day0_atac_NN
#$ -j y
#$ -o ../../../../out/kevin/Writeup6l/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup6l_day0_atac_NN.R