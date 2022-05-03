#!/bin/bash
#$ -N timeAll_saver
#$ -j y
#$ -o ../../../../out/kevin/Writeup4d/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup4d_timeAll_saver.R
