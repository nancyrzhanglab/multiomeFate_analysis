# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5737740/
keygenes_csc <- c("OCT4", "SOX2", "KLF4", 
                  "ALDH1A1", "CD34", "CD47", "CD49f",
                  "CD166", "TNFRSF16", "CD27", "TNFRSF7",
                  "CD105", "CD24", "CD151", "CD38",
                  "CD44", "CD133", "CD166", "CD20", "MS4A1",
                  "CD15", "CD26", "CD90", "CD47",
                  "CD96", "CD29", "CD90", "CD19",
                  "CD117", "CD123", "CEACAM-6", "CD66c", 
                  "CD13", "ALCAM", "TROP1", "GPR49", "CXCR4",
                  "CXCR1", "CX3CR1", "SOX2",
                  "OCT4", "Musashi-1", "Nestin", "BMI1", "CXCL12",
                  "CD45", "NANOG", "TIM3", "REX-1")

# tmp <- sapply(keygenes_csc, function(gene){
#   tmp <- Signac::LookupGeneCoords(
#     object = all_data,
#     gene = gene,
#     assay = "ATAC"
#   )
# })
# table(sapply(tmp, function(i){all(is.null(i))}))