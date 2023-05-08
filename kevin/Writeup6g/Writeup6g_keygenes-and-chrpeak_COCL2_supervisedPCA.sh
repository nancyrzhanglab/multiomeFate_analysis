#!/bin/bash
#$ -N cocl2_spca
#$ -j y
#$ -o ../../../../out/kevin/Writeup6g/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup6g_keygenes-and-chrpeak_COCL2_supervisedPCA.R