#!/bin/bash
#$ -N COCL_ll_b3
#$ -j y
#$ -o ../../../../out/kevin/Writeup6f/qsub/
#$ -l m_mem_free=50G

Rscript --no-save Writeup6f_COCL_loglikelihood_bulk3.R