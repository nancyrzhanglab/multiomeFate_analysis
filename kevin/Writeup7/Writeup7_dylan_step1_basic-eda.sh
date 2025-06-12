#!/bin/bash
#$ -N step1
#$ -j y
#$ -o ../../../../out/kevin/Writeup7/qsub/
#$ -l m_mem_free=50G

Rscript --no-save Writeup7_dylan_step1_basic-eda.R