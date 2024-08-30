#!/bin/bash
#$ -N step7d
#$ -j y
#$ -o ../../../../out/kevin/Writeup10a/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup10a_ppStep7d_chromVar.R
