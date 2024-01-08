#!/bin/bash
#$ -N subset_to
#$ -l m_mem_free=50G
#$ -j y

module load R
Rscript subset_to_DabTram.R