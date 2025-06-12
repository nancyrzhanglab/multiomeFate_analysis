#!/bin/bash
#$ -N lightweight
#$ -j y
#$ -o ../../../../out/kevin/Writeup6o/qsub/
#$ -l m_mem_free=50G

Rscript --no-save Writeup6o_extract_all-data_lightweight.R
