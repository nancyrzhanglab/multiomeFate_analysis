#!/bin/bash
#$ -N icb_fate
#$ -j y
#$ -o ../../../../out/kevin/Writeup3b/qsub
#$ -l m_mem_free=50G

Rscript --no-save w20210731_mouseicb_shell.R
