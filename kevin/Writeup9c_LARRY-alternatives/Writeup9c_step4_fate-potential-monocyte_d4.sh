#!/bin/bash
#$ -N step4_mono
#$ -j y
#$ -o ../../../../out/kevin/Writeup9c/qsub/
#$ -l m_mem_free=10G

Rscript --no-save Writeup9c_step4_fate-potential-monocyte_d4.R
