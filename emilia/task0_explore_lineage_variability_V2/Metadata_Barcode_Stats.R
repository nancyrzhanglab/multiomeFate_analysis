rm(list = ls())
library(ggplot2)
library(gridExtra)
library(tidyverse)

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
            axis.text = element_text(size = 16),
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


# put in the colors
dataset_colors <- c(Absent.day0 = "white",
                    Absent.day10_CIS = "white",
                    Absent.day10_COCL2 = "white",
                    Absent.day10_DABTRAM = "white",
                    Absent.week5_CIS = "white",
                    Absent.week5_COCL2 = "white",
                    Absent.week5_DABTRAM = "white",
                    Present.day0 = "darkgray",
                    Present.day10_CIS = "#FBD08C",
                    Present.day10_COCL2 = "#6DC49C",
                    Present.day10_DABTRAM = "#9D85BE",
                    Present.week5_CIS = "#C96D29",
                    Present.week5_COCL2 = "#0F8241",
                    Present.week5_DABTRAM = "#623594")
dataset_colors <- c(day0 = "darkgray",
                    day10_CIS = "#FBD08C",
                    day10_COCL2 = "#6DC49C",
                    day10_DABTRAM = "#9D85BE",
                    week5_CIS = "#C96D29",
                    week5_COCL2 = "#0F8241",
                    week5_DABTRAM = "#623594")

# =============================================================================
# reading data
# =============================================================================
data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/'
metadat <- read.csv(paste0(data_dir, 'Writeup6b_all-data_metadata.csv'), row.names = 1)

# =============================================================================
# wrangle data
# =============================================================================
metadat$assigned_lineage_presence <- ifelse(is.na(metadat$assigned_lineage), 'Absent', 'Present')
metadat.summary <- metadat %>%
  group_by(dataset, assigned_lineage_presence) %>%
  summarize(n = n())

num.barcode.detected.summary <- metadat %>%
  filter(!is.na(assigned_lineage)) %>%
  group_by(dataset) %>%
  summarize(n = n_distinct(assigned_lineage))
total.barcode <- length(unique(metadat$assigned_lineage))

# =============================================================================
# plot
# =============================================================================

ggplot(metadat.summary, aes(x = dataset, y = n, 
                            fill = interaction(assigned_lineage_presence, dataset))) +
  geom_bar(stat = 'identity', position = 'stack', color = 'black', size = 0.5, width = 0.7) +
  scale_fill_manual(values = dataset_colors) +
  ylim(0, 10500) +
  ylab('Number of cells') +
  xlab('Dataset') +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(metadat.summary, aes(x = dataset, y = n, fill = dataset)) +
  geom_bar(stat = 'identity', position = 'stack', color = 'black', size = 0.5, width = 0.7) +
  scale_fill_manual(values = dataset_colors) +
  ylim(0, 10500) +
  ylab('Number of cells') +
  xlab('Dataset') +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(num.barcode.detected.summary) +
  geom_bar(aes(x = dataset, y = n, fill = dataset), stat = 'identity',color = 'black', size = 0.5, width = 0.7) +
  geom_hline(yintercept = total.barcode, linetype = 'dashed', color = 'black') +
  scale_fill_manual(values = dataset_colors) +
  ylim(0, total.barcode + 100) +
  ylab('Number of barcodes') +
  xlab('Dataset') +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# ===============================================================================================
# ===============================================================================================
# ===============================================================================================

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
out_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task5_FP_for_d10_adapting_cells/'

remove_unassigned_cells <- TRUE

# =============================================================================
# Read data
# =============================================================================
load(paste0(data_dir, 'Writeup10a_data_empty.RData'))

# remove cells with no lineage
if(remove_unassigned_cells) {
  print("Removing cells with no assigned lineage")
  all_data$keep <- !is.na(all_data$assigned_lineage)
  if(any(!all_data$keep)){
    print(paste0("There are ", length(which(!all_data$keep)), " cells being removed"))
    all_data <- subset(all_data, keep == TRUE)
  }
}

metadat <- all_data@meta.data

# =============================================================================
# wrangle data
# =============================================================================
lin.size <- metadat %>%
  group_by(dataset, assigned_lineage) %>%
  summarize(n = n())
lin.size <- spread(lin.size, key = dataset, value = n)

lin.size[is.na(lin.size)] <- 0

lin.size[, c(2:8)] <- log10(lin.size[, c(2:8)] + 1)
lin.size.day10 <- lin.size[, c('assigned_lineage', 'day10_CIS', 'day10_COCL2', 'day10_DABTRAM')]
lin.size.week5 <- lin.size[, c('assigned_lineage', 'week5_CIS', 'week5_COCL2', 'week5_DABTRAM')]


# =============================================================================
# plot
# =============================================================================
p1 <- ggplot(lin.size.day10, aes(x = day10_DABTRAM, y = day10_COCL2)) +
  geom_jitter(width = 0.1, color = '#D3D3D3', shape = 21, fill = "#E0E0E0") +
  stat_cor(size = 7, label.sep = "\n") +
  xlab('DABTRAM') +
  ylab('COCL2') +
  theme_Publication()

p2 <- ggplot(lin.size.day10, aes(x = day10_CIS, y = day10_COCL2)) +
  geom_jitter(width = 0.1, color = '#D3D3D3', shape = 21, fill = "#E0E0E0") +
  stat_cor(size = 7, label.sep = "\n") +
  xlab('CIS') +
  ylab('COCL2') +
  theme_Publication()

p3 <- ggplot(lin.size.day10, aes(x = day10_DABTRAM, y = day10_CIS)) +
  geom_jitter(width = 0.1, color = '#D3D3D3', shape = 21, fill = "#E0E0E0") +
  stat_cor(size = 7, label.sep = "\n") +
  xlab('DABTRAM') +
  ylab('CIS') +
  theme_Publication()


grid.arrange(p1, p2, p3, ncol = 3)


p4 <- ggplot(lin.size.week5, aes(x = week5_DABTRAM, y = week5_COCL2)) +
  geom_jitter(width = 0.1, color = '#D3D3D3', shape = 21, fill = "#E0E0E0") +
  stat_cor(size = 7, label.sep = "\n") +
  xlab('DABTRAM') +
  ylab('COCL2') +
  theme_Publication()

p5 <- ggplot(lin.size.week5, aes(x = week5_CIS, y = week5_COCL2)) +
  geom_jitter(width = 0.1, color = '#D3D3D3', shape = 21, fill = "#E0E0E0") +
  stat_cor(size = 7, label.sep = "\n") +
  xlab('CIS') +
  ylab('COCL2') +
  theme_Publication()

p6 <- ggplot(lin.size.week5, aes(x = week5_DABTRAM, y = week5_CIS)) +
  geom_jitter(width = 0.1, color = '#D3D3D3', shape = 21, fill = "#E0E0E0") +
  stat_cor(size = 7, label.sep = "\n") +
  xlab('DABTRAM') +
  ylab('CIS') +
  theme_Publication()


grid.arrange(p4, p5, p6, ncol = 3)



p7 <- ggplot(lin.size, aes(x = day10_DABTRAM, y = week5_DABTRAM)) +
  geom_jitter(width = 0.1, color = 'black', shape = 21, fill = "#E0E0E0") +
  stat_cor(size = 4, label.sep = "\n", label.x = -0.1, label.y = 2.8) +
  xlab('Day 10') +
  ylab('Week 5') +
  ggtitle('DABTRAM') +
  theme_Publication()


p8 <- ggplot(lin.size, aes(x = day10_COCL2, y = week5_COCL2)) +
  geom_jitter(width = 0.1, color = 'black', shape = 21, fill = "#E0E0E0") +
  stat_cor(size = 4, label.sep = "\n", label.x = -0.1, label.y = 2.8) +
  xlab('Day 10') +
  ylab('Week 5') +
  ggtitle('COCL2') +
  theme_Publication()

p9 <- ggplot(lin.size, aes(x = day10_CIS, y = week5_CIS)) +
  geom_jitter(width = 0.1, color = 'black', shape = 21, fill = "#E0E0E0") +
  stat_cor(size = 4, label.sep = "\n", label.x = -0.1, label.y = 2.8) +
  xlab('Day 10') +
  ylab('Week 5') +
  ggtitle('CIS') +
  theme_Publication()

grid.arrange(p7, p8, p9, ncol = 3)






