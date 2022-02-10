#!/bin/bash
#$ -N t0_fasttopics
#$ -j y
#$ -o ../../../../out/kevin/Writeup4a/qsub/
#$ -l m_mem_free=100G

Rscript --no-save time0_2022_02_09_fasttopics.R
