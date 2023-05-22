#!/bin/bash
#$ -N tcca_RNA-ATAC_DABTRAM
#$ -j y
#$ -o ../../../../out/kevin/Writeup6h/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup6h_tcca_RNA-ATAC_DABTRAM.R