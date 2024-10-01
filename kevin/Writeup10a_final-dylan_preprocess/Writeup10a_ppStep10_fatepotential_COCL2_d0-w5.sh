#!/bin/bash
#$ -N step10_cocl2-d0-w5
#$ -j y
#$ -o ../../../../out/kevin/Writeup10a/qsub/
#$ -l m_mem_free=20G

Rscript --no-save Writeup10a_ppStep10_fatepotential_COCL2_d0-w5.R
