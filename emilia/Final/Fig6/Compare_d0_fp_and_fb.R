rm(list=ls())
library(Seurat)
library(tidyverse)
library(ggpubr)
library(ggdensity)
library(ggplot2)
library(RColorBrewer)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
out_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task5_FP_for_d10_adapting_cells/'
figure_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/figures/MultiomeFate_V3/Fig6/'

remove_unassigned_cells <- TRUE
treatment <- 'DABTRAM'


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
            legend.position = "right",
            legend.box.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "mm"),
            legend.direction = "vertical",
            legend.key.size= unit(0.5, "cm"),
            legend.spacing = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"),
            strip.background=element_rect(colour="#F0F0F0",fill="#F0F0F0"),
            strip.text = element_text(face="bold")
    ))
}

dataset_colors <- c(day0 = "#696969",
                    day10_CIS = "#FBD08C",
                    day10_COCL2 = "#6DC49C",
                    day10_DABTRAM = "#9D85BE",
                    week5_CIS = "#C96D29",
                    week5_COCL2 = "#0F8241",
                    week5_DABTRAM = "#623594")

# ==============================================================================
# Read data
# ==============================================================================
load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_fasttopic_DABTRAM.RData'))
load(paste0(data_dir, 'Writeup10a_data_fatepotential.RData'))

all_data[[paste0("fasttopic.DABTRAM")]] <- all_data_fasttopic_DABTRAM
all_data[[paste0("ft.DABTRAM.umap")]] <- all_data_ft_DABTRAM_umap
# remove cells with no lineage
if(remove_unassigned_cells) {
  print("Removing cells with no assigned lineage")
  all_data$keep <- !is.na(all_data$assigned_lineage)
  if(any(!all_data$keep)){
    print(paste0("There are ", length(which(!all_data$keep)), " cells being removed"))
    all_data <- subset(all_data, keep == TRUE)
  }
}


fp.d0_d10 <- as.data.frame(all_data_fatepotential[[paste0("fatepotential_", treatment, "_d0_d10")]][["cell_imputed_score"]])
fp.d10_w5 <- as.data.frame(all_data_fatepotential[[paste0("fatepotential_", treatment, "_d10_w5")]][["cell_imputed_score"]])

# fate bias from d0 to week5 adapting
df.bias <- read.csv(paste0(out_dir, 'adapting_bias_thres_0_', treatment, '.csv'))

# lineage data
metadat <- all_data@meta.data
metadat.day10_DABTRAM <- metadat[metadat$dataset == 'day10_DABTRAM',]
metadat.week5_DABTRAM <- metadat[metadat$dataset == 'week5_DABTRAM',]

# umap
ft.umap.DABTRAM <- all_data@reductions[[paste0("ft.DABTRAM.umap")]]@cell.embeddings
ft.umap.DABTRAM <- as.data.frame(ft.umap.DABTRAM)
ft.umap.DABTRAM$cell_id <- rownames(ft.umap.DABTRAM)

# ==============================================================================
# Wrangle
# ==============================================================================
metadat <- all_data@meta.data
metadat$cell_id <- rownames(metadat)

colnames(fp.d0_d10) <- c("fatepotential_d0_d10")
fp.d0_d10$cell_id <- rownames(fp.d0_d10)

df <- merge(df.bias, fp.d0_d10, by = 'cell_id')
df <- merge(df, metadat[, c('cell_id', 'assigned_lineage')], by = 'cell_id')

df$imputed_cell_count <- 10**df$fatepotential_d0_d10
df$adapting_cell_count <- 10**df$adaptingFP
df$frac_adapting <- df$adapting_cell_count / df$imputed_cell_count * 100
df$frac_adapting <- ifelse(df$frac_adapting > 100, 100, df$frac_adapting)

df <- df[order(df$frac_adapting, decreasing = F),]
ggplot(df, aes(x = fatepotential_d0_d10, y = bias)) +
  geom_point(aes(fill = adapting_cell_count, size = adapting_cell_count),  shape = 21, color = '#888888', stroke = 0.5) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_hline(yintercept = 0.5, linetype = 'dashed') +
  scale_fill_gradient(low = "#DCDCDC", high = "#7C00FE") +
  # scale_size_continuous(breaks = c(0, 0.1, 0.2, 0.5)) +
  scale_size(range = c(0, 5)) +
  ylim(0, 1) +
  xlim(-1.1, 1) +
  labs(title = treatment,
       y = 'Day0 fate bias to Day10-adapting',
       x = 'Fate potential from Day0 to Day10') +
  theme_Publication()
ggsave(paste0(figure_dir, 'adapting_bias_thres_0_vs_fp_d0_d10_DABTRAM.pdf'), width = 5.2, height = 3.5)


lineage.size <- df %>% 
  group_by(assigned_lineage) %>% 
  summarise(d0.lin.size = n(),
            adaptingCount = sum(adapting_cell_count),
            maxBias = max(bias),
            max.d0_d10_FP = max(fatepotential_d0_d10)) %>% 
  arrange(desc(d0.lin.size))

df.1 <- df[df$assigned_lineage == 'Lin88999', ]
df.1 <- df[df$assigned_lineage == 'Lin111903', ]
df.1 <- df[df$assigned_lineage == 'Lin11293', ]
df.1 <- df[df$assigned_lineage == 'Lin117363', ]
df.1 <- df[df$assigned_lineage == 'Lin105555', ]
ggplot(df.1, aes(x = fatepotential_d0_d10, y = bias)) +
  geom_point(aes(fill = adapting_cell_count, size = adapting_cell_count),  shape = 21, color = '#888888', stroke = 0.5) +
  geom_point(data = df, aes(fill = adapting_cell_count, size = adapting_cell_count),  shape = 21, color = NA, stroke = 0.5) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_hline(yintercept = 0.5, linetype = 'dashed') +
  scale_fill_gradient(low = "#DCDCDC", high = "#7C00FE") +
  # scale_size_continuous(breaks = c(0, 0.1, 0.2, 0.5)) +
  scale_size(range = c(0, 5)) +
  ylim(0, 1) +
  xlim(-1.1, 1) +
  labs(title = paste0(treatment, ': Lin117363'),
       y = 'Day0 fate bias to Day10-adapting',
       x = 'Fate potential from Day0 to Day10') +
  theme_Publication() +
  theme(legend.position = 'none')
# ggsave(paste0(figure_dir, 'adapting_bias_thres_0_vs_fp_d0_d10_DABTRAM_Lin117363.pdf'), width = 3, height = 3)


lineage.size.W5 <- metadat.week5_DABTRAM %>% 
  group_by(assigned_lineage) %>% 
  summarise(W5.lin.size = n()) %>% 
  arrange(desc(W5.lin.size))
lineage.size.D10 <- metadat.day10_DABTRAM %>% 
  group_by(assigned_lineage) %>% 
  summarise(D10.lin.size = n()) %>% 
  arrange(desc(D10.lin.size))

lineage.size <- merge(lineage.size, lineage.size.D10, by = 'assigned_lineage')
lineage.size <- merge(lineage.size, lineage.size.W5, by = 'assigned_lineage')

ft.umap.DABTRAM <- merge(ft.umap.DABTRAM, metadat[, c('cell_id', 'assigned_lineage', 'dataset')], by = 'cell_id')

colnames(fp.d10_w5) <- c('fp.d10_w5')
fp.d10_w5$cell_id <- rownames(fp.d10_w5)

lin.to.plot <- 'Lin79040' # Priming 'Lin79040' 'Lin18546'
lin.to.plot <- 'Lin56667' # Size 'Lin16087'
lin.to.plot <- 'Lin125036' # Plasticity 

df.1 <- df[df$assigned_lineage == lin.to.plot, ]

p1 <- ggplot(df.1, aes(x = fatepotential_d0_d10, y = bias)) +
  geom_point(data = df, color = '#E8E8E8', stroke = 0.5) +
  geom_point(shape = 18, size = 5, color = '#696969', stroke = 0.5) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_hline(yintercept = 0.5, linetype = 'dashed') +
  # scale_fill_gradient(low = "#DCDCDC", high = "#7C00FE") +
  # scale_size_continuous(breaks = c(0, 0.1, 0.2, 0.5)) +
  scale_size(range = c(0, 5)) +
  ylim(0, 1) +
  xlim(-1.1, 1) +
  labs(title = paste0(treatment, ': ', lin.to.plot),
       y = 'Day0 fate bias to Day10-adapting',
       x = 'Fate potential from Day0 to Day10') +
  theme_Publication() +
  theme(legend.position = 'none')

ft.umap.DABTRAM.d10 <-  merge(ft.umap.DABTRAM, fp.d10_w5, by = 'cell_id')
ft.umap.DABTRAM.d10 <- ft.umap.DABTRAM.d10[order(ft.umap.DABTRAM.d10$fp.d10_w5, decreasing = FALSE), ]
ft.umap.DABTRAM.d0.w5 <- ft.umap.DABTRAM[ft.umap.DABTRAM$dataset %in% c('day0', 'week5_DABTRAM'), ]


max_val <- stats::quantile(ft.umap.DABTRAM.d10$fp.d10_w5, 
                           probs = 0.99, 
                           na.rm = TRUE)
min_val <- stats::quantile(ft.umap.DABTRAM.d10$fp.d10_w5, 
                           probs = 0.01, 
                           na.rm = TRUE)
ft.umap.DABTRAM.d10$fp.d10_w5_scaled <- scales::rescale(
  pmax(pmin(ft.umap.DABTRAM.d10$fp.d10_w5, 
            max_val), 
       min_val),
  to = c(-1, 1)
)



p2 <- ggplot(ft.umap.DABTRAM, aes(x = ftDABTRAMumap_1, y = ftDABTRAMumap_2)) +
  geom_point(color = '#E8E8E8', size = 0.5) +
  geom_point(data = ft.umap.DABTRAM.d0.w5[ft.umap.DABTRAM.d0.w5$assigned_lineage == lin.to.plot, ], 
             aes(color = dataset), size = 5, shape = 18) +
  geom_point(data = ft.umap.DABTRAM.d10[ft.umap.DABTRAM.d10$assigned_lineage == lin.to.plot, ],
             aes(fill = fp.d10_w5_scaled), shape = 21, size = 5, stroke = 0) +
  scale_color_manual(values = dataset_colors) +
  scale_fill_gradient2(low = 'blue', high = 'red', mid = 'bisque', midpoint = 0, limits = c(min_val, max_val)) +
  xlab('') +
  ylab('') +
  theme_Publication() +
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = 'none')

p3 <- ggarrange(p1, p2, widths = c(1, 1), ncol = 2)

ggsave(paste0(figure_dir, 'adapting_bias_thres_0_vs_fp_d0_d10_DABTRAM_panel_', lin.to.plot, '.pdf'), p3, width = 5.8, height = 3)


p4 <- ggplot(ft.umap.DABTRAM, aes(x = ftDABTRAMumap_1, y = ftDABTRAMumap_2)) +
  geom_point(color = '#E8E8E8', size = 0.5) +
  geom_point(data = ft.umap.DABTRAM.d0.w5[ft.umap.DABTRAM.d0.w5$assigned_lineage == lin.to.plot, ], 
             aes(color = dataset), size = 5, shape = 18) +
  geom_point(data = ft.umap.DABTRAM.d10[ft.umap.DABTRAM.d10$assigned_lineage == lin.to.plot, ],
             aes(fill = fp.d10_w5_scaled), shape = 21, size = 5, stroke = 0) +
  scale_color_manual(values = dataset_colors) +
  scale_fill_gradient2(low = 'blue', high = 'red', mid = 'bisque', midpoint = 0, limits = c(min_val, max_val)) +
  xlab('') +
  ylab('') +
  theme_Publication() +
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())
ggsave(paste0(figure_dir, 'adapting_bias_thres_0_vs_fp_d0_d10_DABTRAM_panel_', lin.to.plot, '_legend.pdf'), p4, width = 5.8, height = 3)

