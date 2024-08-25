#!/bin/bash
#$ -N step4
#$ -j y
#$ -o ../../../../out/kevin/Writeup10a/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup10a_ppStep4_lineage.R
