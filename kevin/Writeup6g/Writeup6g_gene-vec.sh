#!/bin/bash
#$ -N gene_vec
#$ -j y
#$ -o ../../../../out/kevin/Writeup6g/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup6g_gene-vec.R