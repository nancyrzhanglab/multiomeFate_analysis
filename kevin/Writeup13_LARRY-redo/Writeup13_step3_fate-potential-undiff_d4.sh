#!/bin/bash
#$ -N step3_undiff
#$ -j y
#$ -o ../../../../out/kevin/Writeup13/qsub/
#$ -l m_mem_free=20G

Rscript --no-save Writeup13_step3_fate-potential-undiff_d4.R
