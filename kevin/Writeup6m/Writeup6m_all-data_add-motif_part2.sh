#!/bin/bash
#$ -N all-data_motif_p2
#$ -j y
#$ -o ../../../../out/kevin/Writeup6m/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup6m_all-data_add-motif_part2.R