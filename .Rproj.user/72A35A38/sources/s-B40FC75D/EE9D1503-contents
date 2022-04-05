############################
# 1. QDECR - Generation R 
###########################
# author: Lorenza Dall'Aglio, as based on scripts by Sander Lamballais
# project: brain morphology - ADHD, confounding bias
# dataset: the Generation R Study @9

### Setting Up the Environment ###

# Load Packages

library(QDECR)
library(mice)

# Set paths, dates, open data

## set wd
dir_wd <- "SETYOURPATH"
setwd(dir_wd)

## directories for the QDECR analyses
dir_tmp <- "/dev/shm"
dir_subj <- "SETFREESURFERPATH"

## directory output
dir_out <- "SETYOURPATH"

project_date <- "010321"

## read in the imputed data
dd11 <- readRDS("/PATH/imp30-30_010321.rds")

## set hemispheres to both lh and rh to loop for both later
hemis = c("lh", "rh")


### Select participants with neuro data ###

#loop for IDs to match (ID match bw IDC and ids in the subject MRI folders) 
fsDir <- "SETFREESURFERDIRECTORY"
dd11$data$usable_fs_f9_final <- NA
for (s in 1:length(dd11$data$idc)){
  fPath <- file.path(fsDir, dd11$data$idc[s], 'surf')   #random check for whether ids exist for that neuro. data
  if (file.exists(fPath)){
    dd11$data$usable_fs_f9_final[s] <- 1
  }
}


# check if there are any NAs 
## if there are, you need to select the ids with no NAs only

table(dd11$data$usable_fs_f9_final,useNA="always") #no NAs


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


# Qdecr run 
## here for every model (specified as a formula in the f list), 
# we run a vertex wise analysis, for each hemisphere. 

for (i in seq_along(f)) {
  project <- paste0("M", i, "_", project_date)
  for (hemi in hemis) {
    vw <- qdecr_fastlm(as.formula(f[[i]]), 
                       data = dd11,
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


# Qdecr run 

for (i in seq_along(f)) {
  project <- paste0("M", i, "_", project_date)
  for (hemi in hemis) {
    vw <- qdecr_fastlm(as.formula(f[[i]]), 
                       data = dd11,
                       id = "idc",
                       hemi = hemi,
                       dir_tmp = dir_tmp,
                       dir_out = dir_out,
                       project = project) 
 }
}

