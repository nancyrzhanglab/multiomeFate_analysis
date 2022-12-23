#!/bin/bash
#$ -N coverageplot_track
#$ -j y
#$ -o ../../../../out/kevin/Writeup6b/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup6b_coverageplot_track.R
