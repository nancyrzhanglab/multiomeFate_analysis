#!/bin/bash
#$ -N step10_cis_d0
#$ -j y
#$ -o ../../../../out/kevin/Writeup10a/qsub/
#$ -l m_mem_free=20G

Rscript --no-save Writeup10a_ppStep10_fatepotential_CIS_d0-d10.R
