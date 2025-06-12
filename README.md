The analysis accompanying the codebase https://github.com/nancyrzhanglab/multiomeFate.

To reproduce the analysis, run the code starting at `ppStep1_peakmerging.R`, and proceeding in order of the numbering of `ppStep`. The files and a brief description of each script is as follows:

- `ppStep1_peakmerging.R`: Create one peak set based on all the ATAC 10x runs
- `ppStep2_rna-merge.R`: Create one gene set based on all the RNA 10x runs
- `ppStep3_combine.R`: Combine the RNA and ATAC data across all the 10x runs into one Seurat object
- `ppStep4_lineage.R`: Prepare the lineage count matrix from the 10x runs
- `ppStep5_barcode-assignment.R`: Apply CloneClean to assign cells to a lineage
- `ppStep5b_lineage-plotting.R`: Plot diagnostics based on the lineage assignment
- `ppStep6_qc-step1_compute.R`: Compute QC metrics for the cells
- `ppStep6_qc-step2_threshold.R`: Filter cells based on the QC metrics
- `ppStep7a_fasttopics.R`: Apply FastTopics to create a low-dimensional embedding of the RNA
- `ppStep7b_saver.R`: Apply SAVER to denoise the RNA gene expression
- `ppStep7c_peakvi1_r-to-py.R`: Prepare the data for PeakVI in Python
- `ppStep7c_peakvi2_All.py`, `ppStep7c_peakvi2_CIS.py`, `ppStep7c_peakvi2_COCL2.py`, and `ppStep7c_peakvi2_DABTRAM.py`: Apply PeakVI in Python to create a low-dimensional embedding of the ATAC
- `ppStep7d_chromVar.R`: Apply ChromVAR to create TF motif scores based on the ATAC
- `ppStep8_combining.R`: Combine FastTopics, SAVER, PeakVI, and chromVAR results into the Seurat object
- `ppStep8b_wnn.R`: Apply WNN to create an UMAP combining ATAC and RNA
- `ppStep8c_saver-pca.R`: Apply PCA to the SAVER denoised gene expression
- `ppStep9_basic-plotting.R`: Make initial plots of the data
- `ppStep10_fatepotential_CIS_d0-d10.R`, `ppStep10_fatepotential_CIS_d10-w5.R`, `ppStep10_fatepotential_COCL2_d0-d10.R`, `ppStep10_fatepotential_COCL2_d10-w5.R`, `ppStep10_fatepotential_DABTRAM_d0-d10.R`, and `ppStep10_fatepotential_DABTRAM_d10-w5.R`: Apply CYFER to each pair of consecutive time points for each treatment
- `ppStep11_fatepotential_combine.R`: Put all the CYFER fate potential results into the Seurat object
- `ppStep12_fatepotential_plotting.R`: Make plots based on the CYFER results