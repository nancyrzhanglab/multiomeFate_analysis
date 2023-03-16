#!/bin/bash
#$ -N cis_coverage
#$ -j y
#$ -o ../../../../out/kevin/Writeup6d/qsub/
#$ -l m_mem_free=150G

Rscript --no-save Writeup6d_CIS_chromatin-coverage.R