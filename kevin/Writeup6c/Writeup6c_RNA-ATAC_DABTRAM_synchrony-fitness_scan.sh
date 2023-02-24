#!/bin/bash
#$ -N RNA-ATAC_DABTRAM_synchrony-fitness_scan
#$ -j y
#$ -o ../../../../out/kevin/Writeup6c/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup6c_RNA-ATAC_DABTRAM_synchrony-fitness_scan.R