#!/bin/bash
#$ -N tcca_RNA-geneAct
#$ -j y
#$ -o ../../../../out/kevin/Writeup5a/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup5a_tiltedCCA_tiltedCCA_RNA-geneActivity.R