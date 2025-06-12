#!/bin/bash
#$ -N step3
#$ -j y
#$ -o ../../../../out/kevin/Writeup8/qsub/
#$ -l m_mem_free=50G

Rscript --no-save Writeup8_step3_fasttopics.R
