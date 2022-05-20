#!/bin/bash
#$ -N splicemerge_cocl2
#$ -j y
#$ -o ../../../../out/kevin/Writeup4d/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup4d_timeAll_splicedUnsplice_seuratMerge_COCL2.R
