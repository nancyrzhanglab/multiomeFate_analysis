#!/bin/bash
#$ -N peakvi_all
#$ -j y
#$ -o ../../../../out/kevin/Writeup10a/qsub/
#$ -l m_mem_free=50G

conda activate scvi

conda env list
python --version
conda list

python $HOME/project/Multiome_fate/git/multiomeFate_analysis_kevin/kevin/Writeup10a_final-dylan_preprocess/Writeup10a_ppStep7c_peakvi2_All.py