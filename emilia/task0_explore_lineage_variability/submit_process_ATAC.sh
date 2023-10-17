#!/bin/bash
#$ -N process_ATAC
#$ -l m_mem_free=100G
#$ -j y

module load R
Rscript --vanilla /home/mnt/weka/nzh/team/emiliac/nzhanglab/project/emiliac/lineage_trace/task0_explore_lineage_variability/scripts/process_ATAC.R NT_dICB_D19_InVivo_Rep2

