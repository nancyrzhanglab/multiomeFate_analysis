library(Seurat)
library(ggpubr)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'

theme_Publication<- function(base_size=12, base_family="sans") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            plot.subtitle = element_text(face = "bold", hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold"),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line.x = element_line(colour="black"),
            axis.line.y = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="white"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.box.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "mm"),
            legend.direction = "horizontal",
            legend.key.size= unit(0.5, "cm"),
            legend.spacing = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"),
            strip.background=element_rect(colour="#F0F0F0",fill="#F0F0F0"),
            strip.text = element_text(face="bold")
    ))
}

############################################################################################################
# Read data general
############################################################################################################

load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_fatepotential.RData'))
load(paste0(data_dir, 'Writeup10a_data_saver.RData'))

all_data@misc <- all_data_fatepotential
all_data[["Saver"]] <- all_data_saver
all_data[["Saver.pca"]] <- all_data_saver_pca
all_data[["Saver.umap"]] <- all_data_saver_umap

# ==============================================================================
# Plot correlation
# ==============================================================================

# subset data to day 10
all_data_d10_dabtram <- subset(all_data, dataset == 'day10_DABTRAM')
saver_d10_dabtram <- all_data_d10_dabtram@assays[["Saver"]]@data

saver_d10_dabtram_g <- as.data.frame(saver_d10_dabtram['FOSL1', ])
colnames(saver_d10_dabtram_g) <- c('FOSL1_Saver')

fate_potential_d10_dabtram <- as.data.frame(all_data_fatepotential[["fatepotential_DABTRAM_d10_w5"]][["cell_imputed_score"]])
colnames(fate_potential_d10_dabtram) <- c('fatepotential_DABTRAM_d10_w5')

to_plot <- merge(saver_d10_dabtram_g, fate_potential_d10_dabtram, by = 'row.names')

ggplot(to_plot, aes(x = FOSL1_Saver, y = fatepotential_DABTRAM_d10_w5)) +
  geom_point(color = '#9D85BE') +
  stat_cor(method = "pearson", label.x = 0.1, label.y = 2.5) +
  coord_cartesian(xlim = c(0, 2)) +
  theme_Publication()

# ==============================================================================
# Plot correlation
# ==============================================================================

# subset data to day 10
all_data_d10_cocl2 <- subset(all_data, dataset == 'day10_COCL2')
saver_d10_cocl2 <- all_data_d10_cocl2@assays[["Saver"]]@data

saver_d10_cocl2_g <- as.data.frame(saver_d10_cocl2['TIMP3', ])
colnames(saver_d10_cocl2_g) <- c('TIMP3_Saver')

fate_potential_d10_cocl2 <- as.data.frame(all_data_fatepotential[["fatepotential_COCL2_d10_w5"]][["cell_imputed_score"]])
colnames(fate_potential_d10_cocl2) <- c('fatepotential_COCL2_d10_w5')

to_plot <- merge(saver_d10_cocl2_g, fate_potential_d10_cocl2, by = 'row.names')

ggplot(to_plot, aes(x = TIMP3_Saver, y = fatepotential_COCL2_d10_w5)) +
  geom_point(color = '#6DC49C') +
  stat_cor(method = "pearson", label.x = 0.1, label.y = 2.5) +
  # coord_cartesian(xlim = c(0, 2)) +
  theme_Publication()

# ==============================================================================
# Plot correlation
# ==============================================================================

# subset data to day 10
all_data_d0 <- subset(all_data, dataset == 'day0')
saver_d0<- all_data_d0@assays[["Saver"]]@data

saver_d0_g <- as.data.frame(saver_d0['FOSL1', ])
colnames(saver_d0_g) <- c('FOSL1_Saver')

fate_potential_d0_dabtram <- as.data.frame(all_data_fatepotential[["fatepotential_DABTRAM_d0_d10"]][["cell_imputed_score"]])
colnames(fate_potential_d0_dabtram) <- c('fatepotential_DABTRAM_d0_d10')

to_plot <- merge(saver_d0_g, fate_potential_d0_dabtram, by = 'row.names')
to_plot$fatepotential_DABTRAM_d0_d10 <- as.numeric(as.character(to_plot$fatepotential_DABTRAM_d0_d10))

ggplot(to_plot, aes(x = FOSL1_Saver, y = fatepotential_DABTRAM_d0_d10)) +
  geom_point(color = '#778899') +
  stat_cor(method = "pearson", label.x = 0.1, label.y = 1.2) +
  coord_cartesian(xlim = c(0, 2)) +
  theme_Publication()

