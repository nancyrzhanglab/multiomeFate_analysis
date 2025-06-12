#!/bin/bash
#$ -N step7c
#$ -j y
#$ -o ../../../../out/kevin/Writeup10a/qsub/
#$ -l m_mem_free=50G

Rscript --no-save Writeup10a_ppStep7c_peakvi1_r-to-py.R
