#!/bin/bash
#$ -N barcode_estimation
#$ -j y
#$ -o ../../../../out/kevin/Writeup4f/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup4f_timeAll_barcode_estimation.R
