#!/bin/bash
#$ -N chrom_saver
#$ -j y
#$ -o ../../../../out/kevin/Writeup6b/qsub/
#$ -l m_mem_free=50G
#$ -pe openmp 4

Rscript --no-save Writeup6b_chromatinAct_saver.R
