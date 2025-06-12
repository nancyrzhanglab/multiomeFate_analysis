#!/bin/bash
#$ -N fp_mono-4
#$ -j y
#$ -o ../../../../out/kevin/Writeup11/qsub/
#$ -l m_mem_free=30G

Rscript --no-save Writeup11_step3_fatepotential_monocyte-d4.R
