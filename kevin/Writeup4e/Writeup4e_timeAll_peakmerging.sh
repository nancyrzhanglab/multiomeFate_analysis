#!/bin/bash
#$ -N timeAll_peakmerging
#$ -j y
#$ -o ../../../../out/kevin/Writeup4e/qsub/
#$ -l m_mem_free=150G

Rscript --no-save Writeup4e_timeAll_peakmerging.R
