#!/bin/bash
#$ -N priming_v2
#$ -j y
#$ -o ../../../../out/kevin/Writeup14/qsub/
#$ -l m_mem_free=20G

Rscript --no-save Writeup14_priming_v2.R
