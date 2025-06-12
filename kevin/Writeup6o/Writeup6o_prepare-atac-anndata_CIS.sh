#!/bin/bash
#$ -N anndata_cis
#$ -j y
#$ -o ../../../../out/kevin/Writeup6o/qsub/
#$ -l m_mem_free=50G

Rscript --no-save Writeup6o_prepare-atac-anndata_CIS.R