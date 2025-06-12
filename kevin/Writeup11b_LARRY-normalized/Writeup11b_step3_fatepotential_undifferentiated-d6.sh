#!/bin/bash
#$ -N fp_undiff-6
#$ -j y
#$ -o ../../../../out/kevin/Writeup11b/qsub/
#$ -l m_mem_free=30G

Rscript --no-save Writeup11b_step3_fatepotential_undifferentiated-d6.R
