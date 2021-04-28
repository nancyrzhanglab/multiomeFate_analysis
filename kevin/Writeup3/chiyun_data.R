rm(list=ls())
lis <- readRDS("../../../../../../Downloads/myeloid_traj2_df2000_processed.rds")
names(lis)
head(lis$df_x)
dim(lis$df_x)
head(lis$df_y)
dim(lis$df_y)

mat <- readRDS("../../../../../../Downloads/myeloid_traj2_traj_atac2000_interploation.rds")
dim(mat)
image(.rotate(mat))

list_traj_mat <- readRDS("../../../../../../Downloads/list_traj_mat.rds")
obj_next <- readRDS("../../../../../../Downloads/obj_next.rds")

image(.rotate(list_traj_mat[[1]]))

# start fidgeting around
order_vec <- stats::hclust(stats::dist(t(list_traj_mat[[1]])))$order
image(.rotate(list_traj_mat[[1]][,order_vec]))

##############################
# now for real

rm(list=ls())
tmp <- readRDS("../../../../../../Downloads/myeloid_traj2_df2000_processed.rds")
df_x <- tmp$df_x; df_y <- tmp$df_y
list_traj_mat <- readRDS("../../../../../../Downloads/list_traj_mat.rds")
tmp <- readRDS("../../../../../../Downloads/obj_next.rds")
image(.rotate(tmp$mat_g))

set.seed(10)
obj_next <- prepare_obj_nextcell(df_x, df_y, mat_g = tmp$mat_g, list_traj_mat = list_traj_mat)

set.seed(10)
dat <- generate_data(obj_next)

plot(obj_next$obj_blueprint$vec_colmean)
image(.rotate(dat$obs_x))

############

y_obs <- dat$obs_y[5,]; y_true <- dat$true_y[5,]

image(.rotate(obj_next$ht[["1"]]$list_coef$mat_coef))
