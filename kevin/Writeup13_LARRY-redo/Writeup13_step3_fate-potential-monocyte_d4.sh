#!/bin/bash
#$ -N step3_mono
#$ -j y
#$ -o ../../../../out/kevin/Writeup13/qsub/
#$ -l m_mem_free=20G

Rscript --no-save Writeup13_step3_fate-potential-monocyte_d4.R
