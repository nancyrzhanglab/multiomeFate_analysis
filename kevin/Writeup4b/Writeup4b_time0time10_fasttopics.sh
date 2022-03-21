#!/bin/bash
#$ -N t0t10_fasttopics
#$ -j y
#$ -o ../../../../out/kevin/Writeup4b/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup4b_time0time10_fasttopics.R
