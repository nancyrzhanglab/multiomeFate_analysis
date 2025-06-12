#!/bin/bash
#$ -N step2
#$ -j y
#$ -o ../../../../out/kevin/Writeup13/qsub/
#$ -l m_mem_free=50G

Rscript --no-save Writeup13_step2_fasttopics.R
