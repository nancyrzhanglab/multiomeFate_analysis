#!/bin/bash
#$ -N anndata_lightweight
#$ -j y
#$ -o ../../../../out/kevin/Writeup6o/qsub/
#$ -l m_mem_free=35G

Rscript --no-save Writeup6o_prepare-atac-anndata_lightweight.R