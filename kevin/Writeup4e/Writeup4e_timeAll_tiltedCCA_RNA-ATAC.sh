#!/bin/bash
#$ -N tcca_RNA-ATAC
#$ -j y
#$ -o ../../../../out/kevin/Writeup4e/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup4e_timeAll_tiltedCCA_RNA-ATAC.R