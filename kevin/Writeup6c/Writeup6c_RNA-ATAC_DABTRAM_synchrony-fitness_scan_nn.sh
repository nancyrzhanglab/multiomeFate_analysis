#!/bin/bash
#$ -N RNA-ATAC_DABTRAM_synchrony-fitness_scan_nn
#$ -j y
#$ -o ../../../../out/kevin/Writeup6c/qsub/
#$ -l m_mem_free=50G

Rscript --no-save Writeup6c_RNA-ATAC_DABTRAM_synchrony-fitness_scan_nn.R