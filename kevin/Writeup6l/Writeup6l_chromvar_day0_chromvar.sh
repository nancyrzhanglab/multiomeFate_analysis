#!/bin/bash
#$ -N day0-chromvar
#$ -j y
#$ -o ../../../../out/kevin/Writeup6l/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup6l_chromvar_day0_chromvar.R