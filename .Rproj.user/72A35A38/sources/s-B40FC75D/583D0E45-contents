
#### PACKAGES ####

x <- c("mice", "table1", "tidyverse")
lapply(x, require, character.only = T) 

## setting working directory
setwd("/Users/hannah/Documents/Research/ABCD Research/Confounding project")

import <- readRDS("nda2.0.1.RDS")

# Get var names
varnames <- dput(names(import))
write.csv(varnames, file="varnames.csv")
myvars <- readLines("myvars_toupload.csv")

####################################
# Subset vars from my excel list   #
mydata <- import[myvars]
####################################
# Subset vars to baseline arm 1 and acceptable QC 
mydata <- subset(mydata, mydata$event_name == "baseline_year_1_arm_1" & mydata$fsqc_qc == "accept")
# QC: 462 reject's, 337 NA's
####################################

# adding 'sub-' to ID's
mydata$src_subject_id <- gsub('_', '', mydata$src_subject_id)
mydata$src_subject_id <- paste("sub-", mydata$src_subject_id, sep="")

write.csv(mydata,file="abcd_dataset.csv")

# Table 1
label(mydata$interview_age) <- "Age"
label(mydata$sex_at_birth) <- "Sex"
label(mydata$race_ethnicity) <- "Race/Ethnicity"
label(mydata$household.income) <- "Household Income"
label(mydata$high.educ) <- "Highest Education"
label(mydata$devhx_8_tobacco_p) <- "Tobacco use during preg"
label(mydata$devhx_8_alcohol_p) <- "Alcohol use during preg"
label(mydata$devhx_8_marijuana_p) <- "Marijuana use during preg"
label(mydata$devhx_8_coc_crack_p) <- "Cocaine use during preg"
label(mydata$devhx_8_her_morph_p) <- "Morphine use during preg"
label(mydata$devhx_8_oxycont_p) <- "Oxycontin use during preg"
label(mydata$devhx_3_age_at_birth_mother_p) <- "Mother age at birth"
label(mydata$asr_scr_totprob_t) <- "Primary caregiver psychopathology"
label(mydata$pea_wiscv_tss) <- "Child IQ"

table1(~  + interview_age + sex_at_birth + race_ethnicity  + household.income + high.educ + devhx_8_tobacco_p + devhx_8_alcohol_p + devhx_8_marijuana_p + devhx_8_coc_crack_p + devhx_8_her_morph_p + devhx_8_oxycont_p , data=mydata)

# percent of missingness for each variable
care_vars <- c("interview_age","sex_at_birth" , "race_ethnicity","household.income", "high.educ" , "devhx_8_tobacco_p" , "devhx_8_alcohol_p" , "devhx_8_marijuana_p" , "devhx_8_coc_crack_p", "devhx_8_her_morph_p", "devhx_8_oxycont_p","cbcl_scr_syn_attention_t","devhx_3_age_at_birth_mother_p","asr_scr_totprob_t","pea_wiscv_tss","cbcl_scr_syn_aggressive_t","cbcl_scr_syn_anxdep_t","devhx_12a_born_premature_p","rel_family_id")

x <- apply(mydata, 2, function(x) sum(is.na(x))/length(x)*100)
sort(x)
md.pattern(mydata)

# We run the mice code with 0 iterations 
imp <- mice(mydata, maxit=0)

# Extract predictorMatrix and methods of imputation 
pred = imp$predictorMatrix
meth = imp$method

# Setting values of variables I'd like to leave out to 0 in the predictor matrix

pred_vars <- c("src_subject_id", "site_id_l", "interview_age", "abcd_site", 
               "event_name", "fsqc_qc", "mrif_score", "rel_family_id")
pred[, pred_vars] <- 0



# With this command, we tell mice to impute the data, create 5 datasets, use predM as the predictor matrix and don't print the imputation process. 
imp2 <- mice(mydata, 
             maxit = 5, 
             pred = pred, 
             method = "cart")

# #rename the dataframe to something you like....this is handy when you re-load the mice obj back into R later
final_data <- imp2

# #save out the mice object as an Rdata structure (which can be 'load('fsBully.imp.14sept2018.RData')' later)
saveRDS(final_data, "imputed_dataset.rds",version=2)

xyz <- readRDS("imputed_dataset.rds")
str(xyz)

```

