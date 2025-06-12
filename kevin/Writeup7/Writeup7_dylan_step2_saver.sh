#!/bin/bash
#$ -N step2
#$ -j y
#$ -o ../../../../out/kevin/Writeup7/qsub/
#$ -l m_mem_free=50G
#$ -pe openmp 4

Rscript --no-save Writeup7_dylan_step2_saver.R
