#!/bin/bash
#$ -N day0-atac-subset_lw
#$ -j y
#$ -o ../../../../out/kevin/Writeup6l/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup6l_day0-atac_extract_lightweight.R