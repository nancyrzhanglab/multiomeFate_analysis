#!/bin/bash
#$ -N testrun_d10-w5
#$ -j y
#$ -o ../../../../out/kevin/Writeup6n/qsub/
#$ -l m_mem_free=50G

Rscript --no-save Writeup6n_CIS_testrun_d10-w5.R