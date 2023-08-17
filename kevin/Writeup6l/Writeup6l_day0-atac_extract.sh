#!/bin/bash
#$ -N day0-atac-subset
#$ -j y
#$ -o ../../../../out/kevin/Writeup6l/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup6l_day0-atac_extract.R