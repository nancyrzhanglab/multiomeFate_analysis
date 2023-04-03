#!/bin/bash
#$ -N cis_entropy
#$ -j y
#$ -o ../../../../out/kevin/Writeup6e/qsub/
#$ -l m_mem_free=150G

Rscript --no-save Writeup6e_CIS_peak-entropy.R