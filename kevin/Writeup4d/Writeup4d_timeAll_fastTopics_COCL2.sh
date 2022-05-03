#!/bin/bash
#$ -N timeAll_fastTopics_COCL2
#$ -j y
#$ -o ../../../../out/kevin/Writeup4d/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup4d_timeAll_fastTopics_COCL2.R
