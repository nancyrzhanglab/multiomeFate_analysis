#!/bin/bash
#$ -N chromvar_de_all
#$ -j y
#$ -o ../../../../out/kevin/Writeup6n/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup6n_chromvar_de_all-timepoints.R