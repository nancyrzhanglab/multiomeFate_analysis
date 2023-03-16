#!/bin/bash
#$ -N cocl2_coverage
#$ -j y
#$ -o ../../../../out/kevin/Writeup6d/qsub/
#$ -l m_mem_free=150G

Rscript --no-save Writeup6d_COCL2_chromatin-coverage.R