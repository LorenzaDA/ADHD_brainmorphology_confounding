#### PATHS
dir_wd <- "PATH/TO/YOUR/GENR/WORKING/DIRECTORY"
dir_tmp <- file.path(dir_wd, "dir_tmp")
dir_out <-file.path(dir_wd, "qdecr")

# other paths
FREESURFER_HOME <- Sys.getenv("FREESURFER_HOME")
SUBJECTS_DIR <- file.path(FREESURFER_HOME, "subjects")
FSAVERAGE <- file.path(SUBJECTS_DIR, "fsaverage")

#### OTHER VARS
hemis <- c("lh", "rh")
project_date <- "010321"
hemis <- c("lh","rh")
det <- "sum_att_9m"
out <- c("qdecr_area", "qdecr_volume")

#### INIT
setwd(dir_wd)

#### ANNOTATION
annot_dir <- file.path(FREESURFER_HOME, "subjects", "fsaverage", "label")
annot_files <- paste0(hemis, ".aparc.annot")
annot_paths <- file.path(annot_dir, annot_files)
annot <- lapply(annot_paths, load.annot)
lut <- lapply(annot, `[[`, "LUT")
labs_old <- lapply(annot, `[[`, "vd_label")
labs <- lapply(seq_along(lut), function(i) 
  with(lut[[i]], LUT_labelname[match(labs_old[[i]], LUT_value)]))

#### MODELS

# change according to confounders
conf1 <- c("ethninfv3", "agechildbrainmrif9", "gender")
conf2 <- c(conf1, "educm5", "income5", "age_m_v2")
conf3 <- c(conf2, "msmoke", "can_m", "gsi")
conf4 <- c(conf3, "f0300178")

conf <- list(conf1, conf2, conf3, conf4a)
base_f <- lapply(conf, function(x) paste("~", det, "+", paste(x, collapse = "+")))
f <- lapply(base_f, function(x) paste(out, x))

out2 <- sub("qdecr_", "", out)
base_name <- paste0("_", project_date, ".", out2)

m <- length(f) # number of models
paths_lh <- paths_rh <- list()

for (i in seq_along(out2)) {
  project_names <- expand.grid(".M", seq_len(m), base_name[i]) %>%
    apply(1, paste, collapse = "")
  
  projects_lh <- paste0("lh", project_names)
  projects_rh <- paste0("rh", project_names)
  
  paths_lh[[i]] <- file.path(dir_wd, projects_lh)
  paths_rh[[i]] <- file.path(dir_wd, projects_rh)
}
