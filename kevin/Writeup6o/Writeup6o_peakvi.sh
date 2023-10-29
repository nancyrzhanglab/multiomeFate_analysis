#!/bin/bash
#$ -N peakvi
#$ -j y
#$ -o ../../../../out/kevin/Writeup6o/qsub/
#$ -l m_mem_free=50G

source $HOME/project/Multiome_fate/venv396/bin/activate
python $HOME/project/Multiome_fate/git/multiomeFate_analysis_kevin/kevin/Writeup6o/Writeup6o_peakvi.py