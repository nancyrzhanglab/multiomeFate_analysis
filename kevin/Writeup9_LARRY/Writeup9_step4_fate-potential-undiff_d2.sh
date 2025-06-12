#!/bin/bash
#$ -N step4_undiff_d2
#$ -j y
#$ -o ../../../../out/kevin/Writeup9/qsub/
#$ -l m_mem_free=10G

Rscript --no-save Writeup9_step4_fate-potential-undiff_d2.R
