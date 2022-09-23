#### PACKAGES ####
library(QDECR)
library(mice)

#### INIT ####
# wd
dir_wd <- "PATH/TO/YOUR/WORKING/DIRECTORY"
setwd(dir_wd)

## directories for the QDECR analyses
dir_qdecr <- file.path(working_directory, "qdecr")
if (!dir.exists(dir_qdecr)) dir.create(dir_qdecr)

dir_out <- dir_qdecr
if (!dir.exists(dir_out)) dir.create(dir_out)

dir_data <- file.path(working_directory, "datasets")

dir_tmp <- file.path(dir_wd, "dir_tmp")
dir_subj <- Sys.getenv("SUBJECTS_DIR")
dir_fshome <- Sys.getenv("FREESURFER_HOME")

# set project date
project_date <- "150722"

# set hemispheres to both lh and rh to loop for both later
hemis <- c("lh","rh")

#### MODELS ####

# set determinant, outcome, and covariates
det <- c("cbcl_scr_syn_attention_r", "ksads_14_853_p")
out <- c("qdecr_area", "qdecr_volume", "qdecr_thickness")

# change according to confounders
conf1 <- c("race_ethnicity", "sex_at_birth", "interview_age", "abcd_site")
conf2 <- c(conf1, "high.educ","household.income", "devhx_3_age_at_birth_mother_p")
conf3 <- c(conf2, "devhx_8_tobacco_p", "devhx_8_marijuana_p", "asr_scr_totprob_r")
conf4 <- c(conf3, "pea_wiscv_tss")
conf5 <- c(conf3, "euler_number")
conf <- list(conf1, conf2, conf3, conf4, conf5)
base_f <- lapply(conf, function(x) paste("~", det[1], "+", paste(x, collapse = "+")))
f <- lapply(base_f, function(x) paste(out, x))
            
base_kf <- lapply(conf, function(x) paste("~", det[2], "+", paste(x, collapse = "+")))
kf <- lapply(base_kf, function(x) paste(out, x))
            
#### RUN VERTEX-WISE ANALYSES
            
## read in the imputed data
imp <- readRDS(file.path(dir_data, "imputed_dataset.rds"))
          
# continuous
for (i in seq_along(f)){
  for (j in seq_along(f[[i]])) {
    project <- paste0("M", i, "_", project_date)
    for (hemi in hemis){
      vw <- qdecr_fastlm(as.formula(f[[i]][[j]]),
                         data = imp,
                         id = "src_subject_id",
                         hemi = hemi,
                         dir_subj = dir_subj,
                         dir_tmp = dir_tmp,
                         dir_out = dir_out,
                         project = project,
                         clobber = TRUE)
    }
  }
}

# dichotomous
for (i in seq_along(f)){
  for (j in seq_along(f[[i]])) {
    project <- paste0("M", i, "_ksads_", project_date)
    for (hemi in hemis){
      vw <- qdecr_fastlm(as.formula(kf[[i]][[j]]),
                         data = imp,
                         id = "src_subject_id",
                         hemi = hemi,
                         dir_subj = dir_subj,
                         dir_tmp = dir_tmp,
                         dir_out = dir_out,
                         project = project,
                         clobber = TRUE)
    }
  }
}
        
