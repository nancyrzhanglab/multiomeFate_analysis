#!/bin/bash
#$ -N step7b
#$ -j y
#$ -o ../../../../out/kevin/Writeup10a/qsub/
#$ -l m_mem_free=25G
#$ -pe openmp 4

Rscript --no-save Writeup10a_ppStep7b_saver.R
