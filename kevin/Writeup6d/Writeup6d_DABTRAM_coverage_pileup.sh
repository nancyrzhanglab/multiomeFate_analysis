#!/bin/bash
#$ -N dabtram_coverage-pileup
#$ -j y
#$ -o ../../../../out/kevin/Writeup6d/qsub/
#$ -l m_mem_free=150G

Rscript --no-save Writeup6d_DABTRAM_coverage_pileup.R