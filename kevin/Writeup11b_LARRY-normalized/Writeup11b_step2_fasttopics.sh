#!/bin/bash
#$ -N step2
#$ -j y
#$ -o ../../../../out/kevin/Writeup11b/qsub/
#$ -l m_mem_free=25G

Rscript --no-save Writeup11b_step2_fasttopics.R
