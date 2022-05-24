#!/bin/bash
#$ -N ft_COCL2
#$ -j y
#$ -o ../../../../out/kevin/Writeup4e/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup4e_timeAll_fastTopics_COCL2.R
