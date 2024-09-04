#!/bin/bash
#$ -N fasttopics
#$ -j y
#$ -o ../../../../out/kevin/Writeup11/qsub/
#$ -l m_mem_free=25G

Rscript --no-save Writeup11_step1_fasttopics.R
