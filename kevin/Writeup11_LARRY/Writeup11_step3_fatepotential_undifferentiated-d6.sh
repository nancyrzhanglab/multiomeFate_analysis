#!/bin/bash
#$ -N fp_undif-6
#$ -j y
#$ -o ../../../../out/kevin/Writeup11/qsub/
#$ -l m_mem_free=30G

Rscript --no-save Writeup11_step3_fatepotential_undifferentiated-d6.R
