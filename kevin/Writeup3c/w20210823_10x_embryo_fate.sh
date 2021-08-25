#!/bin/bash
#$ -N embryo_de_fate
#$ -j y
#$ -o ../../../../out/kevin/Writeup3c/qsub
#$ -l m_mem_free=50G

Rscript --no-save w20210823_10x_embryo_fate.R
