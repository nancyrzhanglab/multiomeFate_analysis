#!/bin/bash
#$ -N timeAll_saver
#$ -j y
#$ -o ../../../../out/kevin/Writeup4c/qsub/
#$ -l m_mem_free=25G
#$ -pe openmp 4

Rscript --no-save Writeup4c_timeAll_saver.R
