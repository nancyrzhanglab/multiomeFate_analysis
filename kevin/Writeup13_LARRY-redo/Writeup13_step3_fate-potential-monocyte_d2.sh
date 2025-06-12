#!/bin/bash
#$ -N step3_mono_d2
#$ -j y
#$ -o ../../../../out/kevin/Writeup13/qsub/
#$ -l m_mem_free=20G

Rscript --no-save Writeup13_step3_fate-potential-monocyte_d2.R
