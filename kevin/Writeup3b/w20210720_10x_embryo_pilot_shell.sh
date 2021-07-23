#!/bin/bash
#$ -N embryo_fate
#$ -j y
#$ -o ../../../../out/kevin/Writeup3b/qsub
#$ -l m_mem_free=50G

Rscript --no-save w20210720_10x_embryo_pilot_shell.R
