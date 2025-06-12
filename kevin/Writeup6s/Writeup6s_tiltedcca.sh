#!/bin/bash
#$ -N tiltedcca
#$ -j y
#$ -o ../../../../out/kevin/Writeup6s/qsub/
#$ -l m_mem_free=20G

Rscript --no-save Writeup6s_tiltedcca.R