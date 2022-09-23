#### PACKAGES ####
x <- c("data.table", "foreign", "magrittr", "mice", "QDECR")
lapply(x, require, character.only = T) 
setDTthreads(1)

#### Setting Up the Environment ####

# set project date (for later)
project_date <- "150722"

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


### Euler number ###
aseg_files <- file.path(dir_subj, list.files(dir_subj), "stats", "aseg.stats") 
aseg_files2 <- aseg_files[file.exists(aseg_files)]

euler_number <- sapply(aseg_files2, function(x) fread(cmd = paste("grep SurfaceHoles", x))$V4[3])
euler_idc <- sub(".*(v6.0.0//)(.*)(/stats/aseg.stats)", "\\2", aseg_files2)
euler_data <- data.frame(idc = as.character(euler_idc), euler_number = euler_number)

### Subset ###

imp2 <- lapply(1:imp$m, function(i) {
  complete(imp, i) %>%
    merge(euler_data, all.x = TRUE) %>%
    {.[order(as.numeric(.$idc)), ]} %>%
    subset(!is.na(imp$data$sum_att_9m)) %>%                            # has ADHD data
    subset(mri_consent == "yes") %>%                                   # has MRI consent 
    subset(t1_asset_nii == "NO") %>%                                   # scan type
    subset(usable_fs_f9_final == "useable") %>%                        # FS quality
    subset(qdeclgi_f9 == "has data") %>%                               # QDECR
    subset(braces_mri_f9 == "does not have braces") %>%                # no braces
    subset(exclude_incidental == "include") %>%                        # no incidental findings
    subset(twin == "No") %>%                                           # exclude twins
    subset(!duplicated(mother))                                        # exclude siblings
})

### Run the models in QDECR ###

## Continuous

# MODEL 1: inattention, age, sex, ethnicity
# MODEL 2: maternal education and income and maternal age
# MODEL 3: maternal psychopathology, cannabis use during pregnancy, smoking during pregnancy
# MODEL 4: 3 + child IQ
# MODEL 5: 3 + Euler number

# set determinant, and outcome variable, and the covariates for each model
det <- "sum_att_9m"
outs <- c("qdecr_thickness", "qdecr_area", "qdecr_volume")

# change according to the confounders you want to test and to the N models
conf1 <- c("ethninfv3", "agechildbrainmrif9", "gender")
conf2 <- c(conf1, "educm5", "income5", "age_m_v2")
conf3 <- c(conf2, "msmoke", "can_m", "gsi")
conf4 <- c(conf3, "f0300178")
conf5 <- c(conf3, "euler_number")

# prepare the regression formula
conf <- list(conf1, conf2, conf3, conf4, conf5)
base_f <- lapply(conf, function(x) paste("~", det, "+", paste(x, collapse = " + ")))
f <- lapply(outs, function(out) lapply(base_f, function(x) paste(out, x)))

# QDECR run 
## here for every model (specified as a formula in the f list), 
# we run a vertex wise analysis, for each hemisphere. 
for (i in seq_along(f)) {
  for (j in seq_along(f[[i]])) {
    project <- paste0("M", j, "_genr_", project_date)
    for (hemi in hemis) {
      vw <- qdecr_fastlm(as.formula(f[[i]][[j]]), 
                         data = imp2,
                         id = "idc",
                         hemi = hemi,
                         dir_tmp = dir_tmp,
                         dir_out = dir_out,
                         project = project,
                         clobber = TRUE) 
    }
  }
}
