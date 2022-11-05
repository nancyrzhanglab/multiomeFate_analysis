#!/bin/bash
#$ -N ft_chrAct_CIS
#$ -j y
#$ -o ../../../../out/kevin/Writeup6b/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup6b_chromatinAct_fastTopics_CIS.R