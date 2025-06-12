#!/bin/bash
#$ -N step4_undiff
#$ -j y
#$ -o ../../../../out/kevin/Writeup9d/qsub/
#$ -l m_mem_free=20G

Rscript --no-save Writeup9d_step4_fate-potential-undiff_d4.R
