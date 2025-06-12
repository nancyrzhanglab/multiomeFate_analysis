#!/bin/bash
#$ -N step8
#$ -j y
#$ -o ../../../../out/kevin/Writeup10a/qsub/
#$ -l m_mem_free=50G

Rscript --no-save Writeup10a_ppStep8_combining.R
