#!/bin/bash
#$ -N anndata_cocl2
#$ -j y
#$ -o ../../../../out/kevin/Writeup6o/qsub/
#$ -l m_mem_free=50G

Rscript --no-save Writeup6o_prepare-atac-anndata_COCL2.R