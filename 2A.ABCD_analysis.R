library(QDECR)
library(mice)

# setwd
dir_wd <- "PATH/TO/YOUR/WORKING/DIRECTORY"

setwd(dir_wd)

dir_tmp <- file.path(dir_wd, "dir_tmp")
dir_subj <- Sys.getenv("SUBJECTS_DIR")
dir_out <-file.path(dir_wd, "qdecr")

project_date <- "0522"

# Read data
imp <- readRDS("imputed_dataset.rds")

# set hemispheres to both lh and rh to loop for both later
hemis <- c("lh","rh")

imp$data$src_subject_id_final <- NA
for (s in seq_along(imp$data$src_subject_id)){
  fPath <- file.path(dir_subj, imp$data$src_subject_id[s], "surf")
  if(file.exists(fPath)){
    imp$data$src_subject_id_final[s] <- 1
  }
}

# Subset with no NAs
long1 <- complete(imp, action = "long", include = TRUE)
long2 <- subset(long1, long1$src_subject_id_final == 1)
imp_f <- as.mids(long2)

# QDECR MODELS
n_cores = 10

hemis = c("lh", "rh")

for (hemi in hemis){
  
  a <- qdecr_fastlm(qdecr_area ~ cbcl_scr_syn_attention_r + race_ethnicity  + sex_at_birth + 
                      interview_age + abcd_site + high.educ + household.income + I(devhx_3_age_at_birth_mother_p>29) + 
                      devhx_8_tobacco_p + devhx_8_marijuana_p + I(asr_scr_totprob_r>51),
                    data = imp_f,
                    id ="src_subject_id",
                    hemi = hemi,
                    n_cores = n_cores,
                    dir_tmp = dir_tmp,
                    dir_out = dir_out,
                    project = paste0("bias_", project_date),
                    clobber = TRUE)
  
  a <- qdecr_fastlm(qdecr_volume ~ cbcl_scr_syn_attention_r + race_ethnicity  + sex_at_birth + 
                      interview_age + abcd_site + high.educ + household.income + I(devhx_3_age_at_birth_mother_p>29) + 
                      devhx_8_tobacco_p + devhx_8_marijuana_p + I(asr_scr_totprob_r>51),
                    data = imp_f,
                    id = "src_subject_id",
                    hemi = hemi,
                    n_cores = n_cores,
                    dir_tmp = dir_tmp,
                    dir_out = dir_out,
                    project = paste0("bias_", project_date),
                    clobber = TRUE)
  
}



