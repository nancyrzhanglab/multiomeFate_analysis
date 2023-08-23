#!/bin/bash/
#$ -N homer_cocl2_day10
#$ -j y
#$ -o ../../../../out/kevin/Writeup6l/qsub
#$ -l m_mem_free=100G

cd $HOME/project/Multiome_fate/out/kevin/Writeup6l

module load samtools
module load homer

findMotifsGenome.pl Writeup6l_COCL2-split-by-day10_differential_negpeaks.txt hg38 homer/homer_COCL2-split-by-day10_neg/ -size 200 -mask -bg Writeup6l_COCL2-split-by-day10_differential_bgpeaks.txt
findMotifsGenome.pl Writeup6l_COCL2-split-by-day10_differential_pospeaks.txt hg38 homer/homer_COCL2-split-by-day10_pos/ -size 200 -mask -bg Writeup6l_COCL2-split-by-day10_differential_bgpeaks.txt