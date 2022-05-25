#!/bin/bash
#$ -N treatment_plotting
#$ -j y
#$ -o ../../../../out/kevin/Writeup4e/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup4e_treatment-specific_plotting.R