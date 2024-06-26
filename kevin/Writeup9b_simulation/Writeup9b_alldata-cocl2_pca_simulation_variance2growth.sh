#!/bin/bash
#$ -N var2gro_pca
#$ -j y
#$ -o ../../../../out/kevin/Writeup9b/qsub/
#$ -l m_mem_free=20G

Rscript --no-save Writeup9b_alldata-cocl2_pca_simulation_variance2growth.R
