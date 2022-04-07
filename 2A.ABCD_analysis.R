library(QDECR)
library(mice)

# setwd
dir_wd <- "PATH/TO/YOUR/WORKING/DIRECTORY"

setwd(dir_wd)

dir_tmp <- file.path(dir_wd, "dir_tmp")
dir_subj <- Sys.getenv("SUBJECTS_DIR")
dir_out <-file.path(dir_wd, "qdecr")

project_date <- "0304"

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

#### QDECR MODELS ####

# set determinant, outcome, and covariates
det <- "cbcl_scr_syn_attention_r"
out <- c("qdecr_area", "qdecr_volume")

# change according to confounders
conf1 <- c("race_ethnicity","sex_at_birth","interview_age","abcd_site")
conf2 <- c(conf1,"high.educ","household.income","devhx_3_age_at_birth_mother_p")
conf3 <- c(conf2, "devhx_8_tobacco_p", "devhx_8_marijuana_p", "asr_scr_totprob_r")
conf4 <- c(conf3, "pea_wiscv_tss")
conf <- list(conf1,conf2,conf3,conf4)
base_f <- lapply(conf,function(x)paste("~",det,"+",paste(x,collapse="+")))
f <- lapply(base_f,function(x) paste(out,x))


for (i in seq_along(f)){
  for (j in seq_along(f[[i]])) {
    project <- paste0("M",i,"_",project_date)
    for (hemi in hemis){
      vw <- qdecr_fastlm(as.formula(f[[i]][[j]]),
                         data = imp,
                         id = "src_subject_id",
                         hemi = hemi,
                         dir_tmp = dir_tmp,
                         dir_out = dir_out,
                         project = project,
                         clobber = T,
                         n_cores = 10)
    }
  }
}


