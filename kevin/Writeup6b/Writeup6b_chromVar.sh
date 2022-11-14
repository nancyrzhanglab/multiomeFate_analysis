#!/bin/bash
#$ -N chromVar
#$ -j y
#$ -o ../../../../out/kevin/Writeup6b/qsub/
#$ -l m_mem_free=120G

Rscript --no-save Writeup6b_chromVar.R
