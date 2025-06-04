rm(list = ls())
library(tidyverse)

data_dir <- '~/Downloads/'

# =============================================================================
# Read data
# =============================================================================
data <- read.delim(paste0(data_dir, "2024-10-18 10_51 - registered.csv"), sep=",")
# data <- read.delim(paste0(data_dir, "2024-11-4 - registered.csv"), sep=",")

# =============================================================================
# Wrangle
# =============================================================================
lineage_data <- data %>%
  select(Time, 
         parentId = gfp.blob.Connection.IDs...parentId,
         annotationId = gfp.blob.Connection.IDs...annotationId) %>%
  mutate(across(c(parentId, annotationId), as.character)) %>%
  # filter(Time <= 50) %>%  # Let's start with just first 10 timepoints
  mutate(Time = as.numeric(str_extract(Time, "\\d+"))) %>% 
  arrange(Time)

roots <- lineage_data %>%
  filter(parentId == "-1") %>%
  pull(annotationId)

# group by parentID and concatenate annotationID
division_events <- lineage_data %>%
  filter(parentId != "-1") %>%
  group_by(parentId) %>%
  summarise(annotationId = paste(annotationId, collapse = ","),
            n_child = n()) %>%
  filter(n_child > 1) %>%
  ungroup()


same_cell <- lineage_data %>%
  filter(parentId != "-1") %>%
  group_by(parentId) %>%
  summarise(annotationId = paste(annotationId, collapse = ","),
            n_child = n()) %>%
  filter(n_child == 1) %>%
  ungroup()
same_cell <- as.data.frame(same_cell)
same_cell <- same_cell[, c("parentId", "annotationId")]
rownames(same_cell) <- same_cell$parentId

same_cell_first_id <- setdiff(same_cell$parentId, same_cell$annotationId)
same_cell_list <- list()

same_cell_list <- vector(mode = "list", length = length(same_cell_first_id))
names(same_cell_list) <- same_cell_first_id

for(cell in same_cell_first_id){
  cur_child_cell <- same_cell[same_cell$parentId == cell, "annotationId"]
  same_cell_list[[cell]] <- c(same_cell_list[[cell]], cell)
  while(cur_child_cell %in% same_cell$parentId) {
    same_cell_list[[cell]] <- c(same_cell_list[[cell]], cur_child_cell)
    cur_child_cell <- same_cell[cur_child_cell, "annotationId"]
  }
  same_cell_list[[cell]] <- c(same_cell_list[[cell]], cur_child_cell)
}

same_cell_dict <- vector(mode = "list", length = length(same_cell_list))
names(same_cell_dict) <- names(same_cell_list)
for(cell in names(same_cell_list)){
  same_cell_dict[[cell]] <- tail(same_cell_list[[cell]],1)
}

# division_events$continued_division <- ifelse(length(intersect(str_split_fixed(division_events$annotationId, ',', 2), division_events$parentId)) > 0 , TRUE, FALSE)

division_events$continued_division <- lapply(division_events$annotationId, function(x) {
  length(intersect(str_split_fixed(x, ',', 2), division_events$parentId))
})


# check the alternative IDs
for(i in 1:nrow(division_events)) {
  child_1 <- str_split_fixed(division_events$annotationId[i], ',', 2)[1]
  child_2 <- str_split_fixed(division_events$annotationId[i], ',', 2)[2]
  
  same_cell_list_child_1 <- same_cell_list[[child_1]]
  same_cell_list_child_2 <- same_cell_list[[child_2]]
  
  is_child_1_divided <- FALSE
  child_1_division_cell_id <- NA
  for(c in same_cell_list_child_1) {
    if(c %in% division_events$parentId) {
      is_child_1_divided <- TRUE
      child_1_division_cell_id <- c
    }
  }
  
  is_child_2_divided <- FALSE
  child_2_division_cell_id <- NA
  for(c in same_cell_list_child_2) {
    if(c %in% division_events$parentId) {
      is_child_2_divided <- TRUE
      child_2_division_cell_id <- c
    }
  }
  division_events[i, "child_1"] <- child_1
  division_events[i, "child_2"] <- child_2
  division_events[i, "child_1_divided"] <- paste0(is_child_1_divided, ': ', child_1_division_cell_id)
  division_events[i, "child_2_divided"] <- paste0(is_child_2_divided, ': ', child_2_division_cell_id)
}

# count sub-tree size
count_subtree_size <- function(df, node) {
  # print(node)
  does_node_child1_divide <- df[df$parentId == node, 'child_1_divided']
  does_node_child2_divide <- df[df$parentId == node, 'child_2_divided']
  
  does_node_child1_divide <- as.logical(str_split_fixed(does_node_child1_divide, ': ', 2)[1])
  does_node_child2_divide <- as.logical(str_split_fixed(does_node_child2_divide, ': ', 2)[1])

  child_1_subtree_size <- 1
  child_2_subtree_size <- 1
  
  if (does_node_child1_divide) {
    child_1 <- str_split_fixed(df[df$parentId == node, 'annotationId'], ',', 2)[1]
    child_1_division_cell_id <- str_split_fixed(df[df$parentId == node, 'child_1_divided'], ': ', 2)[2]
    child_1_subtree_size <- count_subtree_size(df, child_1_division_cell_id) + child_1_subtree_size
    
  } 
  
  if(does_node_child2_divide) {
    child_2 <- str_split_fixed(df[df$parentId == node, 'annotationId'], ',', 2)[2]
    child_2_division_cell_id <- str_split_fixed(df[df$parentId == node, 'child_2_divided'], ': ', 2)[2]
    child_2_subtree_size <- count_subtree_size(df, child_2_division_cell_id) + child_2_subtree_size
  }
  return (child_1_subtree_size + child_2_subtree_size)
}

division_events$subtree_size <- sapply(division_events$parentId, function(x) count_subtree_size(division_events, x))


# create a data frame 
division_events.df <- data.frame(
  level1= "C11",
  level2= c( rep("C1020",2), "C283"),
  level3= c( rep("C856", 1), rep("C72", 1), rep("Dead", 1))
)
# transform it to a edge list!
edges_level1_2 <- division_events.df %>% select(level1, level2) %>% unique %>% rename(from=level1, to=level2)
edges_level2_3 <- division_events.df %>% select(level2, level3) %>% unique %>% rename(from=level2, to=level3)
edge_list=rbind(edges_level1_2, edges_level2_3)

# Now we can plot that
mygraph <- graph_from_data_frame( edge_list )
# compute degree for node sizes
V(mygraph)$node_size <- c(2, 2, 2, 8, 46, 10)
V(mygraph)$num_pogeny <- V(mygraph)$node_size
V(mygraph)$node_weight <- c('21', '21', '21', '16', '16', '0')

ggraph(mygraph, layout = 'dendrogram', circular = FALSE) + 
  geom_edge_diagonal() +
  # geom_edge_link(color = "gray") +
  geom_node_point(aes(size = num_pogeny, shape = node_weight), stroke = 2, color = "#0fba56") +
  scale_shape_manual(values = c('21' = 21, '16' = 16, '0' = 4), guide = 'none') +
  scale_size_continuous(range = c(1, 10), breaks = c(2, 10, 20, 40)) +
  geom_node_text(aes(label = name), repel = F, nudge_x = 0.2, size = 4) +
  theme_void() +
  theme(legend.position = 'bottom')
ggsave("~/Downloads/division_events_node11.pdf", width = 4, height = 6)

# create a data frame 
division_events.df2 <- data.frame(
  level1= "C8",
  level2= c( rep("C13",1), rep("C68", 4)),
  level3= c( rep("Dead1", 1), rep("C175", 2), rep("C89", 2)),
  level4= c( rep("Dead2", 1), rep("C290", 1), rep("C346", 1), rep("C482", 1), rep("C159", 1))
)

# transform it to a edge list!
edges_level1_2 <- division_events.df2 %>% select(level1, level2) %>% unique %>% rename(from=level1, to=level2)
edges_level2_3 <- division_events.df2 %>% select(level2, level3) %>% unique %>% rename(from=level2, to=level3)
edges_level3_4 <- division_events.df2 %>% select(level3, level4) %>% unique %>% rename(from=level3, to=level4)
edge_list=rbind(edges_level1_2, edges_level2_3, edges_level3_4)

# Now we can plot that
mygraph <- graph_from_data_frame( edge_list )
# compute degree for node sizes
V(mygraph)$node_size <- c(2, 2, 2, 10, 2, 2, 10, 16, 8, 6, 24)
V(mygraph)$num_pogeny <- V(mygraph)$node_size
V(mygraph)$node_weight <- c('21', '21', '21', '0', '21', '21', '0', '16', '16', '16', '16')

ggraph(mygraph, layout = 'dendrogram', circular = FALSE) + 
  geom_edge_diagonal() +
  # geom_edge_link(color = "gray") +
  geom_node_point(aes(size = num_pogeny, shape = node_weight), stroke = 2, color = "#0fba56") +
  scale_shape_manual(values = c('21' = 21, '16' = 16, '0' = 4), guide = 'none') +
  scale_size_continuous(range = c(1, 8), breaks = c(2, 10, 20, 40)) +
  geom_node_text(aes(label = name), repel = F, nudge_x = 0.4, size = 4) +
  theme_void() +
  theme(legend.position = 'bottom')

ggsave("~/Downloads/division_events_node8.pdf", width = 4, height = 6)
