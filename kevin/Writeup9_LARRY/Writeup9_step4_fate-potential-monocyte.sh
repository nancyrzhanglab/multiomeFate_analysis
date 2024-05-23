#!/bin/bash
#$ -N step4_mono
#$ -j y
#$ -o ../../../../out/kevin/Writeup9/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup9_step4_fate-potential-monocyte.R
