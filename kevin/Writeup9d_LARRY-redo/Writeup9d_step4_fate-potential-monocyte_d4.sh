#!/bin/bash
#$ -N step4_mono
#$ -j y
#$ -o ../../../../out/kevin/Writeup9d/qsub/
#$ -l m_mem_free=20G

Rscript --no-save Writeup9d_step4_fate-potential-monocyte_d4.R
