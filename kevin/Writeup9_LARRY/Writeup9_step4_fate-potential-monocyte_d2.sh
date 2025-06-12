#!/bin/bash
#$ -N step4_mono_d2
#$ -j y
#$ -o ../../../../out/kevin/Writeup9/qsub/
#$ -l m_mem_free=10G

Rscript --no-save Writeup9_step4_fate-potential-monocyte_d2.R
