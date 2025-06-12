#!/bin/bash
#$ -N step2
#$ -j y
#$ -o ../../../../out/kevin/Writeup10a/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup10a_ppStep2_rna-merge.R
