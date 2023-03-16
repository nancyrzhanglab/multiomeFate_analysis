#!/bin/bash
#$ -N dabtram_coverage
#$ -j y
#$ -o ../../../../out/kevin/Writeup6d/qsub/
#$ -l m_mem_free=50G

Rscript --no-save Writeup6d_DABTRAM_chromatin-coverage.R