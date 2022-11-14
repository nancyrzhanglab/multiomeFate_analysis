#!/bin/bash
#$ -N coverageplot
#$ -j y
#$ -o ../../../../out/kevin/Writeup6b/qsub/
#$ -l m_mem_free=50G

Rscript --no-save Writeup6b_coverageplot.R
