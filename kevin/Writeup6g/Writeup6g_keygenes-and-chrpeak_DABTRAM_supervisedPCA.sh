#!/bin/bash
#$ -N dabtram_spca
#$ -j y
#$ -o ../../../../out/kevin/Writeup6g/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup6g_keygenes-and-chrpeak_DABTRAM_supervisedPCA.R