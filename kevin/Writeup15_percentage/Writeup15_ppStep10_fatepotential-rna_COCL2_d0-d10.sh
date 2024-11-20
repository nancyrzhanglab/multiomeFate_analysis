#!/bin/bash
#$ -N step10_cocl2_d0_rna
#$ -j y
#$ -o ../../../../out/kevin/Writeup15/qsub/
#$ -l m_mem_free=20G

Rscript --no-save Writeup15_ppStep10_fatepotential-rna_COCL2_d0-d10.R
