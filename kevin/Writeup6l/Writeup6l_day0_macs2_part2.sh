#!/bin/bash
#$ -N day0-macs2_part2
#$ -j y
#$ -o ../../../../out/kevin/Writeup6l/qsub/
#$ -l m_mem_free=100G

source ../../../../venvMacs2/bin/activate
Rscript --no-save Writeup6l_day0_macs2_part2.R
deactivate