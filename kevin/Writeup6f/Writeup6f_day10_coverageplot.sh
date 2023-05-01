#!/bin/bash
#$ -N day10_coverageplot
#$ -j y
#$ -o ../../../../out/kevin/Writeup6f/qsub/
#$ -l m_mem_free=50G

Rscript --no-save Writeup6f_day10_coverageplot.R
