#!/bin/bash
#$ -N fp_neu-6
#$ -j y
#$ -o ../../../../out/kevin/Writeup11b/qsub/
#$ -l m_mem_free=30G

Rscript --no-save Writeup11b_step3_fatepotential_neutrophil-d6.R
