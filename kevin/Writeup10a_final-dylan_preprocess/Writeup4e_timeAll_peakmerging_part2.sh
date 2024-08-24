#!/bin/bash
#$ -N timeAll_peakmerging
#$ -j y
#$ -o ../../../../out/kevin/Writeup4e/qsub/
#$ -l m_mem_free=300G

Rscript --no-save Writeup4e_timeAll_peakmerging_part2.R
