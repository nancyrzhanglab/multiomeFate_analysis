#!/bin/bash
#$ -N ft_chrAct_all
#$ -j y
#$ -o ../../../../out/kevin/Writeup6b/qsub/
#$ -l m_mem_free=120G

Rscript --no-save Writeup6b_chromatinAct_fastTopics_all.R