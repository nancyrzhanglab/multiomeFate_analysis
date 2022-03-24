#!/bin/bash
#$ -N t0t10_saver
#$ -j y
#$ -o ../../../../out/kevin/Writeup4b/qsub/
#$ -l m_mem_free=25G
#$ -pe openmp 4

Rscript --no-save Writeup4b_time0time10_saver.R
