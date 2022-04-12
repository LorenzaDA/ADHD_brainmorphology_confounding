
#### PACKAGES ####

x <- c("mice", "table1", "tidyverse")
lapply(x, require, character.only = T) 

## setting working directory
setwd("PATH/TO/YOUR/WORKING/DIRECTORY")

import <- readRDS("nda2.0.1.RDS")

# Get var names
write.csv(varnames, file = "varnames.csv")
myvars <- readLines("myvars_toupload.csv")

####################################
# Subset vars from my excel list   #
data <- import[myvars]
####################################
# Subset vars to baseline arm 1 and acceptable QC 
data <- subset(data, 
                 event_name == "baseline_year_1_arm_1" &        # baseline arm 1
                   fsqc_qc == "accept" &                        # acceptable image quality
                   !is.na(cbcl_scr_syn_attention_r) &           # remove missings
                   rel_relationship != "twin" &                 # remove twins
                   rel_relationship != "triplet")               # remove triplets

set.seed(2020)
data <- data[sample(nrow(data)), ] # randomly order data
data <- data[!duplicated(data$rel_family_id), ]

# adding 'sub-' to ID's
data$src_subject_id <- gsub('_', '', data$src_subject_id)
data$src_subject_id <- paste("sub-", data$src_subject_id, sep="")

# remove missing surface data
subid <- read.csv("/PATH/TO/missing.rh.area.qdecr.csv")
data <- left_join(data, subid, by ="src_subject_id")
data <- data %>% filter(yes != 1 | is.na(yes)) %>%
  dplyr::select(-c("yes"))

write.csv(data, file = "abcd_confounding_clean.csv")

# Table 1
label(data$interview_age) <- "Age"
label(data$sex_at_birth) <- "Sex"
label(data$race_ethnicity) <- "Race/Ethnicity"
label(data$household.income) <- "Household Income"
label(data$high.educ) <- "Highest Education"
label(data$devhx_8_tobacco_p) <- "Tobacco use during preg"
label(data$devhx_8_alcohol_p) <- "Alcohol use during preg"
label(data$devhx_8_marijuana_p) <- "Marijuana use during preg"
label(data$devhx_8_coc_crack_p) <- "Cocaine use during preg"
label(data$devhx_8_her_morph_p) <- "Morphine use during preg"
label(data$devhx_8_oxycont_p) <- "Oxycontin use during preg"
label(data$devhx_3_age_at_birth_mother_p) <- "Mother age at birth"
label(data$asr_scr_totprob_r) <- "Primary caregiver psychopathology"
label(data$pea_wiscv_tss) <- "Child IQ"


table1(~ + interview_age + sex_at_birth + race_ethnicity + 
         household.income + high.educ + devhx_8_tobacco_p + devhx_8_marijuana_p, 
       data = data)

## IMPUTATION

# exclude vars
exclVar <- c("src_subject_id","site_id_l","interview_age","abcd_site",
             "event_name",
             "fsqc_qc",
             "mrif_score",
             "rel_family_id",
             "cbcl_scr_syn_attention_r",
             "smri_vol_subcort.aseg_intracranialvolume",
             "rel_relationship")

## Check that the methods is pmm for continuous, logreg for binary, and polyreg for categorical (>2 levels). 

ini <- mice(data, 
            m = 1, 
            pred = quickpred(data, exclude = exclVar, method = "spearman"), 
            seed=2020, 
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
imp <- mice(data, 
            m = nImp, 
            pred = pred, 
            method = meth, 
            seed = 2020, 
            maxit = nIter) 

# #save out the mice object as an Rdata structure (which can be 'load('fsBully.imp.14sept2018.RData')' later)
saveRDS(imp, "imputed_dataset.rds")
