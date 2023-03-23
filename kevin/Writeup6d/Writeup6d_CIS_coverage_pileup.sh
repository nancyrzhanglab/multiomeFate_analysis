#!/bin/bash
#$ -N cis_coverage-pileup
#$ -j y
#$ -o ../../../../out/kevin/Writeup6d/qsub/
#$ -l m_mem_free=150G

Rscript --no-save Writeup6d_CIS_coverage_pileup.R