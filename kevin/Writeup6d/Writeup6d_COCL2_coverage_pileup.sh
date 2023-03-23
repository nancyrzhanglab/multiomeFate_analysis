#!/bin/bash
#$ -N cocl2_coverage-pileup
#$ -j y
#$ -o ../../../../out/kevin/Writeup6d/qsub/
#$ -l m_mem_free=150G

Rscript --no-save Writeup6d_COCL2_coverage_pileup.R