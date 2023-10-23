#!/bin/bash
#$ -N anndata
#$ -j y
#$ -o ../../../../out/kevin/Writeup6o/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup6o_prepare-atac-anndata.R