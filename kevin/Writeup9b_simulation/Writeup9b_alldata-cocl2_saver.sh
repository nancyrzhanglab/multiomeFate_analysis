#!/bin/bash
#$ -N saver
#$ -j y
#$ -o ../../../../out/kevin/Writeup9b/qsub/
#$ -l m_mem_free=25G
#$ -pe openmp 4

Rscript --no-save Writeup9b_alldata-cocl2_saver.R
