#!/bin/bash
#$ -N lsi
#$ -j y
#$ -o ../../../../out/kevin/Writeup6n/qsub/
#$ -l m_mem_free=50G

Rscript --no-save Writeup6n_extract-lsi.R