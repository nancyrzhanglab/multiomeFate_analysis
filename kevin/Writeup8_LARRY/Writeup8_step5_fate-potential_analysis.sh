#!/bin/bash
#$ -N step5
#$ -j y
#$ -o ../../../../out/kevin/Writeup8/qsub/
#$ -l m_mem_free=50G

Rscript --no-save Writeup8_step5_fate-potential_analysis.R
