
#### PACKAGES ####

x <- c("data.table", "mice", "table1", "tidyverse")
lapply(x, require, character.only = TRUE) 

#### INIT ####

# setting working directory
setwd("PATH/TO/YOUR/WORKING/DIRECTORY")

# define SUBJECTS_DIR
dir_subj <- "PATH/TO/SUBJECTS/DIR"

# read in data
import <- readRDS("nda2.0.1.RDS")

# Get var names
write.csv(varnames, file = "varnames.csv")
myvars <- readLines("myvars_toupload.csv")

####################################
# Subset vars from my excel list   #
data <- import[myvars]
####################################

# adding 'sub-' to ID's
data$src_subject_id <- gsub('_', '', data$src_subject_id)
data$src_subject_id <- paste("sub-", data$src_subject_id, sep = "")

# add euler number to data
aseg_files <- file.path(list.files(dir_subj, full.names = TRUE), "stats", "aseg.stats") 

euler_number <- lapply(aseg_files, function(x) {
  suppressWarnings(system(paste("grep 'Measure SurfaceHoles'", x), intern = TRUE))
})
euler_number_fail <- sapply(euler_number, function(x) length(x) == 0)
euler_number <- as.numeric(gsub("[^0-9]", "", euler_number[!euler_number_fail]))
euler_id <- sub(".*(untar//)(.*)(/stats/aseg.stats)", "\\2", aseg_files[!euler_number_fail])

euler_data <- data.frame(src_subject_id = euler_id, euler_number = euler_number)

# test whether surface is missing or not
surf_files <- file.path(list.files(dir_subj, full.names = TRUE), "surf", "rh.area.fwhm10.fsaverage.mgh") 
surf_exists <- file.exists(surf_files)
surf_id <- sub(".*(untar//)(.*)(/surf/rh.area.fwhm10.fsaverage.mgh)", "\\2", surf_files)

surf_data <- data.frame(src_subject_id = surf_id, missing_surf = surf_exists)

data <- merge(as.data.table(data), as.data.table(euler_data), all.x = TRUE)
data <- merge(data, as.data.table(surf_data), all.x = TRUE)

class(data) <- "data.frame"


## SUBSET

if_values <- c("Consider clinical referral", "Consider immediate clinical referral")

data_1 <- data %>%
    subset(event_name == "baseline_year_1_arm_1") %>%                  # baseline arm 1
    subset(!is.na(cbcl_scr_syn_attention_r)) %>%                       # has ADHD data
    subset(fsqc_qc == "accept") %>%                                    # acceptable image quality
    subset(missing_surf == TRUE) %>%                                   # surface
    subset(!mrif_score %in% if_values) %>%                             # no incidental findings
    subset(rel_relationship != "twin") %>%                             # exclude twins
    subset(rel_relationship != "triplet") %>%                          # exclude triplet
    subset(!duplicated(rel_family_id))                                 # exclude siblings


## TABLE

# Table 1
label(data_1$interview_age) <- "Age"
label(data_1$sex_at_birth) <- "Sex"
label(data_1$race_ethnicity) <- "Race/Ethnicity"
label(data_1$household.income) <- "Household Income"
label(data_1$high.educ) <- "Highest Education"
label(data_1$devhx_8_tobacco_p) <- "Tobacco use during preg"
label(data_1$devhx_8_alcohol_p) <- "Alcohol use during preg"
label(data_1$devhx_3_age_at_birth_mother_p) <- "Mother age at birth"
label(data_1$asr_scr_totprob_r) <- "Primary caregiver psychopathology"
label(data_1$pea_wiscv_tss) <- "Child IQ"

table1::table1(~ + interview_age + sex_at_birth + race_ethnicity + 
         household.income + high.educ + devhx_8_tobacco_p + devhx_8_marijuana_p, 
       data = data_1)

# exclude vars
exclVar <- c("src_subject_id", "site_id_l", "interview_age", "abcd_site",
             "event_name",
             "fsqc_qc",
             "mrif_score",
             "rel_family_id",
             "cbcl_scr_syn_attention_r",
             "smri_vol_subcort.aseg_intracranialvolume",
             "rel_relationship", "euler_number")

## Check that the methods is pmm for continuous, logreg for binary, and polyreg for categorical (>2 levels). 

ini <- mice(data_1, 
            m = 1, 
            seed = 2020, 
            pred = quickpred(data, exclude = exclVar, method = "spearman"), 
            maxit = 0)

pred <- ini$pred
meth <- ini$meth

for (i in exclVar){
  meth[i] <- ""
  pred[i,] <- 0
  pred[,i] <- 0
}

## Set parameters for the imputation
nImp <- 5
nIter <- 10

### Save files
imp <- mice(data_1, 
            m = nImp, 
            pred = pred, 
            method = meth, 
            seed = 2020, 
            maxit = nIter) 
           
# #save out the mice object as an Rdata structure
saveRDS(imp, "datasets/imputed_dataset.rds")
