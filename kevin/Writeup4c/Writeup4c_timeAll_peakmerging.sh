#!/bin/bash
#$ -N timeAll_peakmerging
#$ -j y
#$ -o ../../../../out/kevin/Writeup4c/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup4c_timeAll_peakmerging.R
