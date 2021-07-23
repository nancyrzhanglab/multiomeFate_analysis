#!/bin/bash
#$ -N embryo_fate
#$ -j y
#$ -o ../../../../out/Writeup14b/qsub/
#$ -l m_mem_free=50G

Rscript --no-save 20210720_10x_embryo_pilot_shell.R
