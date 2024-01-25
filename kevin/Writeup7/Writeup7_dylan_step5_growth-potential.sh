#!/bin/bash
#$ -N step5
#$ -j y
#$ -o ../../../../out/kevin/Writeup7/qsub/
#$ -l m_mem_free=50G

Rscript --no-save Writeup7_dylan_step5_growth-potential.R
