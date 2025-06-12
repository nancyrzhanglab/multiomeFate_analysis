#!/bin/bash
#$ -N plastic
#$ -j y
#$ -o ../../../../out/kevin/Writeup14/qsub/
#$ -l m_mem_free=20G

Rscript --no-save Writeup14_plastic.R
