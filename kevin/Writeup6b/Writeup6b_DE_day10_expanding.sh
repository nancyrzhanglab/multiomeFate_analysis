#!/bin/bash
#$ -N de_day10
#$ -j y
#$ -o ../../../../out/kevin/Writeup6b/qsub/
#$ -l m_mem_free=50G

Rscript --no-save Writeup6b_DE_day10_expanding.R
