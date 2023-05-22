#!/bin/bash
#$ -N gene-peak-matching
#$ -j y
#$ -o ../../../../out/kevin/Writeup6h/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup6h_gene-peak-matching.R