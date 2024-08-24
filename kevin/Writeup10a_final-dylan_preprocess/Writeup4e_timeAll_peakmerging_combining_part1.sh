#!/bin/bash
#$ -N combining_part1
#$ -j y
#$ -o ../../../../out/kevin/Writeup4e/qsub/
#$ -l m_mem_free=300G

Rscript --no-save Writeup4e_timeAll_peakmerging_combining_part1.R
