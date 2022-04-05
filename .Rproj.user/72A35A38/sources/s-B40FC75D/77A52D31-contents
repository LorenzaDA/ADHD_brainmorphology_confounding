#### PACKAGES ####
library(QDECR)
library(mice)

#### Setting Up the Environment ####

# set project date (for later)
project_date <- "010321"

## set wd
working_directory <- "PATH/TO/YOUR/WORKING/DIRECTORY"

## directories for the QDECR analyses
dir_qdecr <- file.path(working_directory, "qdecr")
if (!dir.exists(dir_qdecr)) dir.create(dir_qdecr)

dir_out <- dir_qdecr
if (!dir.exists(dir_out)) dir.create(dir_out)

dir_data <- file.path(working_directory, "datasets")

dir_tmp <- "/dev/shm"
dir_subj <- Sys.getenv("SUBJECTS_DIR")
dir_fshome <- Sys.getenv("FREESURFER_HOME")

## read in the imputed data
imp <- readRDS(file.path(dir_data, paste0("imp30-30_", project_date, ".rds")))

## set hemispheres to both lh and rh to loop for both later
hemis <- c("lh", "rh")

### Run the models in QDECR ###

## Continuous

# MODEL 1: inattention, age, sex, ethnicity
# MODEL 2: maternal education and income and maternal age
# MODEL 3: maternal psychopathology, cannabis use during pregnancy, smoking during pregnancy
# MODEL 4: child IQ

# set determinant, and outcome variable, and the covariates for each model
## Area 
det <- "sum_att_9m"
out <- "qdecr_area"

# change according to the confounders you want to test and to the N models
conf1 <- c("ethninfv3", "agechildbrainmrif9", "gender")
conf2 <- c(conf1, "educm5", "income5", "age_m_v2")
conf3 <- c(conf2, "msmoke", "can_m", "gsi")
conf4 <- c(conf3, "f0300178")

# prepare the regression formula
conf <- list(conf1, conf2, conf3, conf4)
base_f <- lapply(conf, function(x) paste("~", det, "+", paste(x, collapse = " + ")))
f <- lapply(base_f, function(x) paste(out, x))

# QDECR run 
## here for every model (specified as a formula in the f list), 
# we run a vertex wise analysis, for each hemisphere. 
for (i in seq_along(f)) {
  project <- paste0("M", i, "_", project_date)
  for (hemi in hemis) {
    vw <- qdecr_fastlm(as.formula(f[[i]]), 
                       data = imp,
                       id = "idc",
                       hemi = hemi,
                       dir_tmp = dir_tmp,
                       dir_out = dir_out,
                       project = project) 
  }
  
}

## Volume

# set outcome to volume now 
out <- "qdecr_volume"

# adapt formulas with the new out 
f <- lapply(base_f, function(x) paste(out, x))

# QDECR run 
for (i in seq_along(f)) {
  project <- paste0("M", i, "_", project_date)
  for (hemi in hemis) {
    vw <- qdecr_fastlm(as.formula(f[[i]]), 
                       data = imp,
                       id = "idc",
                       hemi = hemi,
                       dir_tmp = dir_tmp,
                       dir_out = dir_out,
                       project = project) 
  }
}