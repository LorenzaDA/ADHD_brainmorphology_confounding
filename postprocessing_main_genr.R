#### PACKAGES
library(dplyr)
library(ggplot2)
library(magrittr)
library(mice)
library(QDECR)
library(stringr)

#### Font for (some of the) plots
font_fam <- "serif"

#### Rerun Freeview (yes/no) [set to FALSE if you just want to rerun the quick parts of the code]
rerun_freeview <- TRUE

#### PRE INIT
main_dir <- "PATH/TO/YOUR/GENR/WORKING/DIRECTORY"

#### FUNCTIONS
source("postprocessing_functions.R")

#### STUDY PARAMETERS
source("postprocessing_parameters_genr.R")

#### LOOP 
combs <- combn(m, 2, simplify = FALSE)
fx <- list()

# structure: outcome -> hemisphere -> model
sem <- betam <- ocn <- list(list(list(), list()), list(list(), list()))

# structure: outcome -> model
roi_m <- list(list(), list())
roi_med <- list(list(), list())
roi_s <- list(list(), list())

# structure: outcome -> wide/long
roi_list <- list(list(), list())

# limits for region figures
lims <- list(
  list(c(-0.005, 0.001), c(-40, 40), c(-40, 40)),
  list(c(-0.015, 0.003), c(-40, 40), c(-40, 40))
)

# load ocn + beta maps  
for (qq in seq_along(out2)) {
  for (i in seq_along(hemis)) {
    for (j in seq_len(m)) {
      to_load <- c(paths_lh[[qq]][j], paths_rh[[qq]][j])[i]
      vw <- qdecr_load(to_load)
      vw <- qdecr_update_path(vw, FREESURFER_HOME, SUBJECTS_DIR)
      ocn[[qq]][[i]][[j]] <- qdecr_read_ocn_mask(vw, det)
      betam[[qq]][[i]][[j]] <- qdecr_read_coef(vw, det)$x
      sem[[qq]][[i]][[j]] <- qdecr_read_se(vw, det)$x
    }
  }
}

# create sign cluster .mgh files with overlay (up to model j)
for (qq in seq_along(out2)) {
  for (i in seq_along(hemis)) {
    ocn_temp <- ocn[[qq]][[i]]
    ocn_all <- rep(0, length(ocn_temp[[1]]))
    for (j in seq_len(m)) {
      ocn_all[ocn_temp[[j]]] <- j
      path_temp <- file.path(dir_out, paste0(hemis[i], base_name[qq], "_figure_3_overlay_to_model", j, ".mgh"))
      save.mgh(as_mgh(ocn_all), path_temp)
    }
  }
}

# create data frames for sign cluster plots
summary.table <- lapply(seq_along(out2), function(qq) {
  temp <- data.frame(model = paste0("M", 1:m),
                     n_vertices_sum = c(sapply(ocn[[qq]][[1]], function(x) sum(x > 0)),
                                        sapply(ocn[[qq]][[2]], function(x) sum(x > 0))),
                     hemisphere = rep(hemis, each = m)
  )
  temp$model2 <- dplyr::recode(temp$model, M1 = "1",
                               M2 = "2",
                               M3 = "3",
                               M4 = "3 \n+ IQ")
  temp
})

# create sign cluster .mgh files with overlay from model 3 to other subsequent models
for (qq in seq_along(out2)) {
  for (i in seq_along(hemis)) {
    ocn_temp <- ocn[[qq]][[i]]
    ocn_all <- rep(0, length(ocn_temp[[1]]))
    ocn_all[ocn_temp[[3]]] <- 1

    ocn_all2 <- ocn_all
    ocn_all2[ocn_temp[[4]] & ocn_temp[[3]]] <- 2
    ocn_all2[ocn_temp[[4]] & !ocn_temp[[3]]] <- 3
    path_temp <- file.path(dir_out, paste0(hemis[i], base_name[qq], "_special_model_3_overlay_to_model4.mgh"))
    save.mgh(as_mgh(ocn_all2), path_temp)
  }
}


# create delta maps for each model comparison
for (qq in seq_along(out2)) {
  for (i in seq_along(hemis)) {
    for (j in combs) {
      map1 <- betam[[qq]][[i]][[j[1]]]
      map2 <- betam[[qq]][[i]][[j[2]]]
      map_x <- (map2 / map1 - 1) * 100 # percentage change from one map to the other
      file_temp <- paste0(hemis[i], base_name[[qq]], "_delta_map_M", j[1], "_M", j[2], ".mgh")
      path_temp <- file.path(dir_out, file_temp)
      save.mgh(as_mgh(map_x), path_temp)
    }
  }
}

# create region-wise data
for (qq in seq_along(out2)) {
  
  # get betas per region
  for (j in seq_len(m)) {
    roi_m[[j]] <- aggregate(c(betam[[qq]][[1]][[j]], betam[[qq]][[2]][[j]]), by = list(c(labs[[1]], labs[[2]])), mean, na.rm = TRUE)[1:35, 2]
    roi_med[[j]] <- aggregate(c(betam[[qq]][[1]][[j]], betam[[qq]][[2]][[j]]), by = list(c(labs[[1]], labs[[2]])), median, na.rm = TRUE)[1:35, 2]
    roi_s[[j]] <- aggregate(c(betam[[qq]][[1]][[j]], betam[[qq]][[2]][[j]]), by = list(c(labs[[1]], labs[[2]])), sd, na.rm = TRUE)[1:35, 2]
  }
  roi_n <- aggregate(c(betam[[qq]][[1]][[1]], betam[[qq]][[2]][[1]]), by = list(c(labs[[1]], labs[[2]])), mean)[, 1]
  
  # wide format
  roi_df <- data.frame(names = roi_n,
                       model_1 = roi_m[[1]],
                       model_2 = roi_m[[2]],
                       model_3 = roi_m[[3]],
                       s_1 = roi_s[[1]],
                       s_2 = roi_s[[2]],
                       s_3 = roi_s[[3]])
  temp_roi_names <- names(roi_df)[-1]
  roi_df$model_12 <- roi_df$model_2 - roi_df$model_1
  roi_df$model_13 <- roi_df$model_3 - roi_df$model_1
  roi_df$model_23 <- roi_df$model_3 - roi_df$model_2
  roi_df$deltaperc12 <- roi_df$model_12 / roi_df$model_1 * 100
  roi_df$deltaperc13 <- roi_df$model_13 / roi_df$model_1 * 100
  roi_df$deltaperc23 <- roi_df$model_23 / roi_df$model_2 * 100
  roi_hilo <- order(roi_df$model_1)
  
  # long format
  roi_df2 <- rbind(data.frame(names = roi_n, model = roi_m[[1]], s = roi_s[[1]], time = 1),
                   data.frame(names = roi_n, model = roi_m[[2]], s = roi_s[[2]], time = 2),
                   data.frame(names = roi_n, model = roi_m[[3]], s = roi_s[[3]], time = 3)
  )
  roi_df2$names <- as.factor(roi_df2$names)
  roi_df2$names <- factor(roi_df2$names, levels = levels(roi_df2$names)[roi_hilo])
  roi_df2$time <- as.factor(roi_df2$time)
  roi_df2$time <- factor(roi_df2$time, levels = rev(levels(roi_df2$time)))
  roi_df2$ylo <- roi_df2$model - roi_df2$s
  roi_df2$yhi <- roi_df2$model + roi_df2$s
  
  # format for correlations
  roi_cor <- data.frame(names = roi_n,
                        mean1 = roi_m[[1]],
                        mean2 = roi_m[[2]],
                        mean3 = roi_m[[3]],
                        mean4 = roi_m[[4]],
                        mean5 = roi_m[[5]],
                        mean6 = roi_m[[6]],
                        median1 = roi_med[[1]], 
                        median2 = roi_med[[2]], 
                        median3 = roi_med[[3]], 
                        median4 = roi_med[[4]], 
                        median5 = roi_med[[5]], 
                        median6 = roi_med[[6]]
  )
  
  roi_list[[qq]] <- list(roi_df, roi_df2, roi_cor)
}

#### FIGURES

## Figure 1
# file names + paths for overlay plots we want
overl_name <- expand.grid(hemis, "_", project_date, ".",c("area", "volume"), "_figure_3_overlay_to_model1.mgh" ) %>%
  apply(1, paste, collapse = "")
overl_paths <- file.path(dir_out, overl_name)

# some settings
overlay_threshold <- c(0.1, 2.5)
zoom <- 1
cs_size <- 16
head_size <- 70
text_size <- 70

# generate tiff files for left and right hemisphere
if (rerun_freeview) {
  quick_snaps(overl_paths[1], hemis[1], dir_out, SUBJECTS_DIR, overlay_threshold, overlay_method = "linearopaque")
  quick_snaps(overl_paths[2], hemis[2], dir_out, SUBJECTS_DIR, overlay_threshold, overlay_method = "linearopaque")
  quick_snaps(overl_paths[3], hemis[1], dir_out, SUBJECTS_DIR, overlay_threshold, overlay_method = "linearopaque")
  quick_snaps(overl_paths[4], hemis[2], dir_out, SUBJECTS_DIR, overlay_threshold, overlay_method = "linearopaque")
}

# create composite images for left and right, then combine
ofig1 <- quick_compose(overl_paths[1], hemis[1], cs = FALSE, text = "GENR", text_size = text_size)
ofig2 <- quick_compose(overl_paths[2], hemis[2], cs = FALSE)
ofig3 <- quick_compose(overl_paths[3], hemis[1], cs = FALSE, text = "GENR", text_size = text_size)
ofig4 <- quick_compose(overl_paths[4], hemis[2], cs = FALSE)

fx$fig_1 <- list()
fx$fig_1$row2 <- magick::image_append(c(ofig1, ofig2))
fx$fig_1$row4 <- magick::image_append(c(ofig3, ofig4))

## Figure 3

# file names + paths for overlay plots we want
overl_name <- expand.grid(hemis, "_", project_date, ".",c("area", "volume"), "_figure_3_overlay_to_model3.mgh" ) %>%
  apply(1, paste, collapse = "")
overl_paths <- file.path(dir_out, overl_name)

# some settings
overlay_threshold <- c(0.1, 2.5)
zoom <- 1
cs_size <- 16
head_size <- 70
text_size <- 70

# generate tiff files for left and right hemisphere
if (rerun_freeview) {
  quick_snaps(overl_paths[1], hemis[1], dir_out, SUBJECTS_DIR, overlay_threshold, overlay_method = "linearopaque")
  quick_snaps(overl_paths[2], hemis[2], dir_out, SUBJECTS_DIR, overlay_threshold, overlay_method = "linearopaque")
  quick_snaps(overl_paths[3], hemis[1], dir_out, SUBJECTS_DIR, overlay_threshold, overlay_method = "linearopaque")
  quick_snaps(overl_paths[4], hemis[2], dir_out, SUBJECTS_DIR, overlay_threshold, overlay_method = "linearopaque")
}

# create composite images for left and right, then combine
ofig1 <- quick_compose(overl_paths[1], hemis[1], cs = FALSE, text = "GENR", text_size = text_size)
ofig2 <- quick_compose(overl_paths[2], hemis[2], cs = FALSE)
ofig3 <- quick_compose(overl_paths[3], hemis[1], cs = FALSE, text = "GENR", text_size = text_size)
ofig4 <- quick_compose(overl_paths[4], hemis[2], cs = FALSE)

fx$fig_3 <- list()
fx$fig_3$row2 <- magick::image_append(c(ofig1, ofig2))
fx$fig_3$row4 <- magick::image_append(c(ofig3, ofig4))

## FIGURE 4 (outcome -> hemisphere -> model)

mask_col <- lapply(c("lh", "rh"), function(x) {
  system.file("extdata", paste0(x, ".fsaverage.cortex.mask.mgh"), package = "QDECR") %>%
    load.mgh %>%
    `$`(x)
}) %>%
  unlist %>%
  as.logical

data4 <- data.frame(
  vd = c(betam[[1]][[1]][[1]], betam[[1]][[2]][[1]], betam[[1]][[1]][[2]], betam[[1]][[2]][[2]], betam[[1]][[1]][[3]], betam[[1]][[2]][[3]],
         betam[[2]][[1]][[1]], betam[[2]][[2]][[1]], betam[[2]][[1]][[2]], betam[[2]][[2]][[2]], betam[[2]][[1]][[3]], betam[[2]][[2]][[3]]),
  vs = c("no", "yes")[c(ocn[[1]][[1]][[1]], ocn[[1]][[2]][[1]], ocn[[1]][[1]][[1]], ocn[[1]][[2]][[1]], ocn[[1]][[1]][[1]], ocn[[1]][[2]][[1]],
                        ocn[[2]][[1]][[1]], ocn[[2]][[2]][[1]], ocn[[2]][[1]][[1]], ocn[[2]][[2]][[1]], ocn[[2]][[1]][[1]], ocn[[2]][[2]][[1]]) + 1],
  mask_col = mask_col,
  hemi = factor(rep(hemis, each = length(betam[[1]][[1]][[1]]), times = 6)),
  model = factor(rep(1:3, each = length(betam[[1]][[1]][[1]]), times = 4)),
  modality = factor(rep(out2, each = length(betam[[1]][[1]][[1]]) * 6))
)

data4_se <- data.frame(
  vd = c(sem[[1]][[1]][[1]], sem[[1]][[2]][[1]], sem[[1]][[1]][[2]], sem[[1]][[2]][[2]], sem[[1]][[1]][[3]], sem[[1]][[2]][[3]],
         sem[[2]][[1]][[1]], sem[[2]][[2]][[1]], sem[[2]][[1]][[2]], sem[[2]][[2]][[2]], sem[[2]][[1]][[3]], sem[[2]][[2]][[3]]),
  vs = c("no", "yes")[c(ocn[[1]][[1]][[1]], ocn[[1]][[2]][[1]], ocn[[1]][[1]][[1]], ocn[[1]][[2]][[1]], ocn[[1]][[1]][[1]], ocn[[1]][[2]][[1]],
                        ocn[[2]][[1]][[1]], ocn[[2]][[2]][[1]], ocn[[2]][[1]][[1]], ocn[[2]][[2]][[1]], ocn[[2]][[1]][[1]], ocn[[2]][[2]][[1]]) + 1],
  mask_col = mask_col,
  hemi = factor(rep(hemis, each = length(sem[[1]][[1]][[1]]), times = 6)),
  model = factor(rep(1:3, each = length(sem[[1]][[1]][[1]]), times = 4)),
  modality = factor(rep(out2, each = length(sem[[1]][[1]][[1]]) * 6))
)

fx$fig_4$b <- ggplot(subset(data4, mask_col & vd != 0 & modality == "area"), aes(x = vd, group = interaction(model, vs), color = model, linetype = vs)) + 
  geom_density(n = 500, size = 1.2) + 
  theme_classic() + 
  ylab("Density") + 
  xlab("Coefficient") +
  xlim(-0.006, 0.001) + 
  #guides(color=guide_legend("Significant in M1")) +
  theme(text = element_text(size = 12, family = font_fam),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  theme(legend.position = "none")

fx$fig_4$d <- ggplot(subset(data4, mask_col & vd != 0 & modality == "volume"), aes(x = vd, group = interaction(model, vs), color = model, linetype = vs)) + 
  geom_density(n = 500, size = 1.2) + 
  theme_classic() + 
  ylab("Density") + 
  xlab("Coefficient") +
  xlim(-0.03, 0.01) + 
  #guides(color=guide_legend("Significant in M1")) +
  theme(text = element_text(size = 12, family = font_fam),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  theme(legend.position = "none")

fx$fig_4$b_se <- ggplot(subset(data4_se, mask_col & vd != 0 & modality == "area"), aes(x = vd, group = interaction(model, vs), color = model, linetype = vs)) + 
  geom_density(n = 500, size = 1.2) + 
  theme_classic() + 
  ylab("Density") + 
  xlab("Coefficient") +
  #guides(color=guide_legend("Significant in M1")) +
  theme(text = element_text(size = 12, family = font_fam),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  theme(legend.position = "none")

fx$fig_4$d_se <- ggplot(subset(data4_se, mask_col & vd != 0 & modality == "volume"), aes(x = vd, group = interaction(model, vs), color = model, linetype = vs)) + 
  geom_density(n = 500, size = 1.2) + 
  theme_classic() + 
  ylab("Density") + 
  xlab("Coefficient") +
  #guides(color=guide_legend("Significant in M1")) +
  theme(text = element_text(size = 12, family = font_fam),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  theme(legend.position = "none")

## FIGURE 6

# file names + paths for overlay plots we want
overl_name <- expand.grid(hemis, "_", project_date, ".",c("area", "volume"), "_special_model_3_overlay_to_model4.mgh" ) %>%
  apply(1, paste, collapse = "")
overl_paths <- file.path(dir_out, overl_name)

# some settings
overlay_threshold <- c(0.1, 2.5)
zoom <- 1
cs_size <- 16
head_size <- 70
text_size <- 70

# generate tiff files for left and right hemisphere
if (rerun_freeview) {
  quick_snaps(overl_paths[1], hemis[1], dir_out, SUBJECTS_DIR, overlay_threshold, overlay_method = "linearopaque")
  quick_snaps(overl_paths[2], hemis[2], dir_out, SUBJECTS_DIR, overlay_threshold, overlay_method = "linearopaque")
  quick_snaps(overl_paths[3], hemis[1], dir_out, SUBJECTS_DIR, overlay_threshold, overlay_method = "linearopaque")
  quick_snaps(overl_paths[4], hemis[2], dir_out, SUBJECTS_DIR, overlay_threshold, overlay_method = "linearopaque")
}

# create composite images for left and right, then combine
ofig1 <- quick_compose(overl_paths[1], hemis[1], cs = FALSE, text = "GENR", text_size = text_size)
ofig2 <- quick_compose(overl_paths[2], hemis[2], cs = FALSE)
ofig3 <- quick_compose(overl_paths[3], hemis[1], cs = FALSE, text = "GENR", text_size = text_size)
ofig4 <- quick_compose(overl_paths[4], hemis[2], cs = FALSE)

fx$fig_6 <- list()
fx$fig_6$row2 <- magick::image_append(c(ofig1, ofig2))
fx$fig_6$row4 <- magick::image_append(c(ofig3, ofig4))


## Supplementary Figure 2

# file names + paths for delta plots we want
delta_name <- expand.grid(hemis, "_", project_date, ".",c("area", "volume"), "_delta_map_M1_M3.mgh" ) %>%
  apply(1, paste, collapse = "")
delta_paths <- file.path(dir_out, delta_name)

# some settings
overlay_threshold <- c(0.01, 30)
overlay_method <- "linear"
zoom <- 1
cs_size <- 20
head_size <- 70
text_size <- 70
cs_value <- as.character(c(paste0("+", overlay_threshold[2]), -overlay_threshold[2]))
cs_value <- paste0(cs_value, "%")

# generate tiff files for left and right hemisphere
if (rerun_freeview) {
  quick_snaps(delta_paths[1], hemis[1], dir_out, SUBJECTS_DIR, overlay_threshold, inverted = TRUE, zoom = zoom, overlay_method = overlay_method)
  quick_snaps(delta_paths[2], hemis[2], dir_out, SUBJECTS_DIR, overlay_threshold, inverted = TRUE, zoom = zoom, overlay_method = overlay_method)
  quick_snaps(delta_paths[3], hemis[1], dir_out, SUBJECTS_DIR, overlay_threshold, inverted = TRUE, zoom = zoom, overlay_method = overlay_method)
  quick_snaps(delta_paths[4], hemis[2], dir_out, SUBJECTS_DIR, overlay_threshold, inverted = TRUE, zoom = zoom, overlay_method = overlay_method)
}

# create composite images for left and right, then combine
dfig1 <- quick_compose(delta_paths[1], hemis[1], cs = FALSE, text = "GENR", text_size = text_size)
dfig2 <- quick_compose(delta_paths[2], hemis[2], cs = TRUE, cs_size = cs_size, cs_value = cs_value, cs_offset = c("+5+18", "+8+238"))
dfig3 <- quick_compose(delta_paths[3], hemis[1], cs = FALSE, text = "GENR", text_size = text_size)
dfig4 <- quick_compose(delta_paths[4], hemis[2], cs = TRUE, cs_size = cs_size, cs_value = cs_value, cs_offset = c("+5+18", "+8+238"))

fx$fig_s2 <- list()
fx$fig_s2$row2 <- magick::image_append(c(dfig1, dfig2))
fx$fig_s2$row4 <- magick::image_append(c(dfig3, dfig4))

#### OUT

out <- list(betam = betam, ocn = ocn, roi_list = roi_list, fx = fx)
out_path <- file.path(main_dir, "genr.rds")
saveRDS(out, out_path)

fx_genr <- fx

