#!/bin/bash
#$ -N var2gro_saver
#$ -j y
#$ -o ../../../../out/kevin/Writeup9b/qsub/
#$ -l m_mem_free=30G

Rscript --no-save Writeup9b_alldata-cocl2_saver_simulation_variance2growth.R
