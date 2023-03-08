#!/bin/bash
#$ -N peak_counter_onpeak
#$ -j y
#$ -o ../../../../out/kevin/Writeup6c/qsub/
#$ -l m_mem_free=50G

Rscript --no-save Writeup6c_peak_counter_onpeak.R