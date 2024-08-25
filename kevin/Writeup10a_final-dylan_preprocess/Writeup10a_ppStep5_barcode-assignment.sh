#!/bin/bash
#$ -N step5
#$ -j y
#$ -o ../../../../out/kevin/Writeup10a/qsub/
#$ -l m_mem_free=30G

Rscript --no-save Writeup10a_ppStep5_barcode-assignment.R
