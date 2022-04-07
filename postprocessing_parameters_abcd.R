#### PATHS
dir_wd <- "PATH/TO/YOUR/ABCD/WORKING/DIRECTORY"
dir_tmp <- file.path(dir_wd, "dir_tmp")
dir_out <-file.path(dir_wd, "qdecr")

# other paths
FREESURFER_HOME <- Sys.getenv("FREESURFER_HOME")
SUBJECTS_DIR <- file.path(FREESURFER_HOME, "subjects")
FSAVERAGE <- file.path(SUBJECTS_DIR, "fsaverage")

setwd(dir_tmp)

############## ANNOTATION
hemis <- c("lh", "rh")
annot_dir <- file.path(FREESURFER_HOME, "subjects", "fsaverage", "label")
annot_files <- paste0(hemis, ".aparc.annot")
annot_paths <- file.path(annot_dir, annot_files)

annot <- lapply(annot_paths, load.annot)

lut <- lapply(annot, `[[`, "LUT")
labs_old <- lapply(annot, `[[`, "vd_label")

labs <- lapply(seq_along(lut), function(i) 
  with(lut[[i]], LUT_labelname[match(labs_old[[i]], LUT_value)]))

#####################

project_date <- "0304"
hemis <- c("lh","rh")
det <- "cbcl_scr_syn_attention_r"
out <- c("qdecr_area", "qdecr_volume")

### QDECR MODELS ####

# change according to confounders
conf1 <- c("race_ethnicity","sex_at_birth","interview_age")
conf2 <- c(conf1,"high.educ","household.income","devhx_3_age_at_birth_mother_p")
conf3 <- c(conf2, "devhx_8_tobacco_p", "devhx_8_marijuana_p", "asr_scr_totprob_r")
conf5 <- c(conf3, "pea_wiscv_tss")
conf <- list(conf1, conf2, conf3, conf4)
base_f <- lapply(conf, function(x) paste("~", det, "+", paste(x, collapse = "+")))
f <- lapply(base_f, function(x) paste(out, x))

out2 <- sub("qdecr_", "", out)
base_name <- paste0("_", project_date, ".", out2)

m <- length(f) # number of models
paths_lh <- paths_rh <- list()

for (i in seq_along(out2)) {
  project_names <- expand.grid(".M", seq_len(m), base_name[i]) %>%
    apply(1, paste, collapse = "")
  
  project_names[m] <- paste0(".M", m, ".", out2[i])
  
  projects_lh <- paste0("lh", project_names)
  projects_rh <- paste0("rh", project_names)
  
  mid_dir <- if (out2[i] == "area") "dir_out_area" else "dir_out_volume"
  mid_dir <- rep(mid_dir, m)
  mid_dir[m] <- "dir_out"
  
  paths_lh[[i]] <- file.path(dir_wd, mid_dir, projects_lh)
  paths_rh[[i]] <- file.path(dir_wd, mid_dir, projects_rh)
}