#!/bin/bash
#$ -N gene_vec
#$ -j y
#$ -o ../../../../out/kevin/Writeup6f/qsub/
#$ -l m_mem_free=40G

Rscript --no-save Writeup6f_gene-vec.R