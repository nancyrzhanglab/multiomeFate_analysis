#!/bin/bash
#$ -N testrun_d0-d10
#$ -j y
#$ -o ../../../../out/kevin/Writeup6n/qsub/
#$ -l m_mem_free=50G

Rscript --no-save Writeup6n_CIS_testrun_d0-d10.R