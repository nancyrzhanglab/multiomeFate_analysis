#!/bin/bash
#$ -N step7a
#$ -j y
#$ -o ../../../../out/kevin/Writeup10a/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup10a_ppStep7a_fasttopics.R
