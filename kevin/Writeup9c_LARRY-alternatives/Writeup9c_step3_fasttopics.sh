#!/bin/bash
#$ -N step3
#$ -j y
#$ -o ../../../../out/kevin/Writeup9c/qsub/
#$ -l m_mem_free=50G

Rscript --no-save Writeup9c_step3_fasttopics.R
