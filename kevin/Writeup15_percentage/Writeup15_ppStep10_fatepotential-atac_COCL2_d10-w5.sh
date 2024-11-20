#!/bin/bash
#$ -N step10_cocl2_d10_atac
#$ -j y
#$ -o ../../../../out/kevin/Writeup15/qsub/
#$ -l m_mem_free=20G

Rscript --no-save Writeup15_ppStep10_fatepotential-atac_COCL2_d10-w5.R
