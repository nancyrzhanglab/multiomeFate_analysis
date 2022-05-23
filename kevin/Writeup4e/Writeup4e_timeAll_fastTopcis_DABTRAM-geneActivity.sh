#!/bin/bash
#$ -N fastTopics_DAB-activity
#$ -j y
#$ -o ../../../../out/kevin/Writeup4e/qsub/
#$ -l m_mem_free=150G

Rscript --no-save Writeup4e_timeAll_fastTopcis_DABTRAM-geneActivity.R
