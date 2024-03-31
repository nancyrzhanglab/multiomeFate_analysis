#!/bin/bash
#$ -N step4b
#$ -j y
#$ -o ../../../../out/kevin/Writeup8/qsub/
#$ -l m_mem_free=50G

Rscript --no-save Writeup8_step4b_fate-potential.R
