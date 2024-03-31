#!/bin/bash
#$ -N step5b
#$ -j y
#$ -o ../../../../out/kevin/Writeup8/qsub/
#$ -l m_mem_free=50G

Rscript --no-save Writeup8_step5b_fate-potential_analysis.R
