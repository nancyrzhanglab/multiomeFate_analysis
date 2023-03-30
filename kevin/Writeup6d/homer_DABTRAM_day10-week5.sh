#!/bin/bash/
#$ -N homer_dt_d10-w5
#$ -j y
#$ -o ../../../../out/kevin/Writeup6d/qsub
#$ -l m_mem_free=100G

cd $HOME/project/Multiome_fate/out/kevin/Writeup6d

module load samtools
module load homer

findMotifsGenome.pl Writeup6d_DABTRAM_day10-week5_differential_negpeaks.txt hg38 homer/homer_DABTRAM_day10-week5_neg/ -size 200 -mask -bg Writeup6d_DABTRAM_day10-week5_differential_bgpeaks.txt
findMotifsGenome.pl Writeup6d_DABTRAM_day10-week5_differential_pospeaks.txt hg38 homer/homer_DABTRAM_day10-week5_pos/ -size 200 -mask -bg Writeup6d_DABTRAM_day10-week5_differential_bgpeaks.txt