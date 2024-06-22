#!/bin/bash
#$ -N step4_neutro
#$ -j y
#$ -o ../../../../out/kevin/Writeup9c/qsub/
#$ -l m_mem_free=10G

Rscript --no-save Writeup9c_step4_fate-potential-neutrophil_d4.R
