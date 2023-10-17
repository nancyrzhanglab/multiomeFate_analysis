#!/bin/bash
#$ -N subset_data
#$ -l m_mem_free=100G
#$ -j y

module load gcc
module load R
Rscript subset_to_invivo_data.R 