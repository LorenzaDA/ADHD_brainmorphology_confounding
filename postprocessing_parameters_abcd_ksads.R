#### PATHS
dir_wd <- "PATH/TO/YOUR/KSADS/WORKING/DIRECTORY"
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

project_date <- "0208"
hemis <- c("lh", "rh")
det <- "ksads_14_853_p1"
out <- c("qdecr_area", "qdecr_volume")

### QDECR MODELS ####

# change according to confounders
conf1 <- c("race_ethnicity","sex_at_birth","interview_age", "abcd_site")
conf2 <- c(conf1,"high.educ","household.income","devhx_3_age_at_birth_mother_p")
conf3 <- c(conf2, "devhx_8_tobacco_p", "devhx_8_marijuana_p", "asr_scr_totprob_r")
conf4 <- c(conf3, "pea_wiscv_tss")
conf <- list(conf1, conf2, conf3, conf4)
base_f <- lapply(conf, function(x) paste("~", det, "+", paste(x, collapse = "+")))
f <- lapply(base_f, function(x) paste(out,x))

out2 <- sub("qdecr_", "", out)
base_name <- paste0("_", project_date, ".", out2)

m <- length(f) # number of models
paths_lh <- paths_rh <- list()

for (i in seq_along(out2)) {
  project_names <- expand.grid(".M", seq_len(m), base_name[i]) %>%
    apply(1, paste, collapse = "")
  
  projects_lh <- paste0("lh", project_names)
  projects_rh <- paste0("rh", project_names)
  
  mid_dir <- rep("dir_out_2022", m)
  
  paths_lh[[i]] <- file.path(dir_wd, mid_dir, projects_lh)
  paths_rh[[i]] <- file.path(dir_wd, mid_dir, projects_rh)
}