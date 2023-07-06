#!/bin/bash
#$ -N DABTRAM_entropy_tfidf
#$ -j y
#$ -o ../../../../out/kevin/Writeup6j/qsub/
#$ -l m_mem_free=50G

Rscript --no-save Writeup6j_DABTRAM_day0_entropy_tfidf-normalization.R