#!/bin/bash
#$ -N gene-peak-matching
#$ -j y
#$ -o ../../../../out/kevin/Writeup6m/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup6m_gene-peak-matching_tss.R
