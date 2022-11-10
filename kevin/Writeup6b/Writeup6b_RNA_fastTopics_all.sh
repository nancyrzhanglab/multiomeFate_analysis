#!/bin/bash
#$ -N ft_rna_all
#$ -j y
#$ -o ../../../../out/kevin/Writeup6b/qsub/
#$ -l m_mem_free=120G

Rscript --no-save Writeup6b_RNA_fastTopics_all.R