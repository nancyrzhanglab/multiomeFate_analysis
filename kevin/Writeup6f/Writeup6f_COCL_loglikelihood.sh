#!/bin/bash
#$ -N COCL_loglikelihood
#$ -j y
#$ -o ../../../../out/kevin/Writeup6f/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup6f_COCL_loglikelihood.R