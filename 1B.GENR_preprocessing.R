
#### FUNCTIONS ####

#### facred: Factor level reduction; merges factor levels (mer) into a new level (new)
## Variables: 
#     x: A factor.
#     mer: Two or more levels that need to be merged.
#     new: The new level name. If NULL, the new name will be `paste(mer, collapse = clps)`.
#     clps: The `collapse` argument for `new`. This is only relevant if `new` is not null.

facred <- function(x, mer, new = NULL, clps = " || "){
  le <- levels(x)
  print(le)
  le2 <- `names<-`(as.list(le), le)
  ind <- sapply(le2, function(y) y %in% mer) %>% unname
  if (is.null(new)) new <- paste(mer, collapse = clps) %>% as.character
  names(le2)[ind] <- new
  levels(x) <- le2
  x
}

#### PACKAGES ####
x <- c("data.table", "foreign", "magrittr", "mice")
lapply(x, require, character.only = T) 

###################
# 0. Data preparation

# Project: brain morphology - ADHD, confounding
# Dataset: the Generation R Study @9
# Authors: Lorenza Dall'Aglio, as based on scripts by Sander Lamballais-Tessenshon
# NB previous steps: 
# data missings (missings (i.e. 888,999)) and missing system have been transformed into missing system in SPSS for readability in R


###################
# STEP 1: Read in and merge the data

# set directories - not specified here for privacy 

working_directory <- "PATH/TO/YOUR/WORKING/DIRECTORY"
dir_data <- file.path(working_directory, "datasets")
dir_descriptives <- file.path(working_directory, "descriptives")

if (!dir.exists(dir_data)) dir.create(dir_data)
if (!dir.exists(dir_descriptives)) dir.create(dir_descriptives)

# get the names of all datasets in the folder 
files_idm <- c("CHILD-ALLGENERALDATA_29012018.sav",
               "130815_ GEDRAGSGROEP Maternal smoking pregnancy.sav",
               "GR1003-BSI D1_22112016.sav")
files_idc <- c("Cannabis_variable.sav", 
               "CHILDCBCL9_10082016.sav", 
               "CHILDMRI_Core_data_distribution_28112017.sav",
               "CHILDSONIQ5_06042012.sav", 
               "F9_MRI_Freesurfer_global_2016_12_09.sav")

files_genr <- c(files_idm, files_idc)
paths_genr <- file.path(dir_data, files_genr)

# read in all files into a list
data <- lapply(paths_genr, function(x)
  x %>%
    foreign::read.spss(to.data.frame = TRUE) %>%
    data.table::as.data.table() %>%
    {data.table::setnames(., tolower(names(.)))}
)

# order by idm and idc, respectively
data[seq_along(files_idm)] <- 
  lapply(data[seq_along(files_idm)], data.table::setorder, idm)
data[-seq_along(files_idm)] <- 
  lapply(data[-seq_along(files_idm)], data.table::setorder, idc)

####
# STEP 2: Merge data
####

# Merge general data, smoking during pregnancy, drug use, BSI by idm (id of the mother) 
# (these datasets do not contain information on the IDC - id of the child, contrarily to other datasets)

# merge
data <- Reduce(function(d1, d2) merge(d1, d2, by = intersect(names(d1), names(d2)), all.x = TRUE), data) 

####
# STEP 3: Select the variables we need (also for imputation, to improve it)
####

# variables needed for inclusion/exclusion criteria: 

# consent, T1 asset nii, braces, incidental findings, freesurfer useability var, qdecLgi, twin, mother (pregnancy).
# variables for the models: 
# inattention, age, sex, ethnicity
# maternal education and income and maternal age
# maternal psychopathology, cannabis use during pregnancy, smoking during pregnancy
# child IQ

# NB also add any variables which might be useful for the imputation. 
# For example, we considered maternal education at birth 
# because this might help us fill in the missingness for maternal education at age 5

vars <- c("idc", "idm", "sum_att_9m", "qdeclgi_f9", "t1_asset_nii", 
          "usable_fs_f9_final", "braces_mri_f9", "mother", "exclude_incidental",
          "mri_consent", "twin", "ethninfv3", "gender", "agechildbrainmrif9", 
          "age_m_v2", "educm5", "income5", "gestbir", "weight",
          "msmoke", "can_m", "gsi", "f0300178", "educm","knicr_tbv_f9", "etiv_f9")

data <- data[, vars, with = FALSE] 

data[, idc := as.character(idc)]
data[, idm := as.character(idm)]
data[, mother := as.character(mother)]

####
# STEP 4: Flow chart
#### 

# INCLUSION: Select individuals with inattention data & sufficient quality MRI data

data <- data %>%
  subset(!is.na(sum_att_9m)) %>%                                     # has ADHD data
  subset(mri_consent == "yes") %>%                                   # has MRI consent 
  subset(t1_asset_nii == "NO") %>%                                   # scan type
  subset(usable_fs_f9_final == "useable") %>%                        # FS quality
  subset(qdeclgi_f9 == "has data") %>%                               # QDECR
  subset(braces_mri_f9 == "does not have braces") %>%                # no braces
  subset(exclude_incidental == "include") %>%                        # no incidental findings
  subset(twin == "No") %>%                                           # exclude twins
  subset(!duplicated(mother))                                        # exclude siblings

class(data) <- "data.frame"

## Final sample = 2531 children 

####
# STEP 5: Collapse levels in categorical variables with many categories
####

# Check the initial levels of the variables

levels(data$ethninfv3)
levels(data$educm5)
levels(data$income5)
levels(data$educm)
levels(data$can_m)

# Specify how you want to recode the variables

ulist <- list(
  list(c("ethninfv3"), #ethnicity into western and non western
       c("African","American, non western","Asian, non western","Cape Verdian","Dutch Antilles",
         "Indonesian","Moroccan","Surinamese-Creole", "Surinamese", "Surinamese-Hindustani", "Surinamese-Unspecified","Turkish"),
       "Non-western"),
  list(c("ethninfv3"),
       c('Dutch',"American,western","Asian, western","European","Oceanie"),
       "Western"),
  
  list("educm5", #education (measure at age 5) into low, medium and high
       c("no education finished", "primary", "secondary, phase 1"),
       "Low"),
  list("educm5",
       "secondary, phase 2",
       "Medium"),
  list("educm5",
       c("higher, phase 1", "higher, phase 2"),
       "High"),
  
  list("educm", #education (measure at birth) into low,medium, high
       c("no education finished", "primary", "secondary, phase 1"), 
       "Low"), 
  list("educm",
       "secondary, phase 2",
       "Medium"),
  list("educm",
       c("higher, phase 1", "higher, phase 2"),
       "High"),
  
  list("income5", # income into low, medium, high (based on Dutch standards) 
       c("Less than € 800", "€ 800-1200", "€ 1200-1600", "€ 1600-2000"),
       "Low"), 
  list("income5", 
       c("€ 2000-2400", "€ 2400-2800", "€ 2800-3200"),
       "Medium"), 
  list("income5", 
       c("€ 3200-4000", "€ 4000-4800", "€ 4800-5600", "More than € 5600"), 
       "High"),
  
  list("can_m", # cannabis into no exposure or exposure during pregnancy
       # (not considering lifetime exposure to mimic ABCD)
       c("No exposure", "Cannabis exposure before pregnancy"),
       "No"), 
  list("can_m", 
       c("Cannabis exposure during pregnancy"),
       "Yes"))

# recode the variables based on the specification above
for (l1 in seq_len(length(ulist))) {
  temp <- ulist[[l1]]
  for (l2 in seq_len(length(temp[[1]]))){
    out <- temp[[1]][l2]
    if (out %in% names(data)){
      old <- temp[[2]]
      new <- temp[[3]]
      data[,out] <- facred(data[,out], old, new)
    }
  }
}

####
# STEP 7: Prepare data for imputation
####

# Get the N & % missingness

missings_n <- apply(data, 2, function(col) {sum(is.na(col))})

write.csv(missings_n,
          file.path(dir_descriptives, "missings_N.csv"), 
          row.names = F)

missings_percent <- apply(data, 2, function(col) {((sum(is.na(col)))/length(col)*100)})

write.csv(missings_percent,
          file.path(dir_descriptives, "missing_stats_percent.csv"), 
          row.names = F)

# Get summary stats

summarystats <- summary(data) 

write.csv(summarystats,
          file.path(dir_descriptives, "descriptive_stats.csv"), 
          row.names = F)


# Get the sd for Suppl T 2
gsisd <- sd(data$gsi, na.rm = T)
age_m_v2_sd <- sd(data$age_m_v2, na.rm = T)
gestbir_sd <- sd(data$gestbir, na.rm = T) 
iq_sd <- sd(data$f0300178, na.rm = T)
age_sd <- sd(data$agechildbrainmrif9, na.rm = T)
sd_all <- cbind(gsisd, age_m_v2_sd, gestbir_sd, iq_sd, age_sd)  

write.csv(sd_all,
          file.path(dir_descriptives, "sds.csv"), 
          row.names = F)


####
# STEP 8: impute
####

### Variables 

# specify the variables that can aid the prediction (i.e. help predict the missingness in other variables)
predictors_for_imputation <- c("age_m_v2", "educm", "income5", "msmoke", "can_m", "sum_att_9m", "knicr_tbv_f9")

# specify the variables that you want to impute. It's important to add in the predictors for imputation 
## this is important because if our predictors have NAs, our output of the imputation will also present NAs 
impvars <- c("gestbir", "weight", "educm5", "f0300178", "gsi", "ethninfv3", "income5", "msmoke", predictors_for_imputation)  

#dryrun mice (page 35 mice guide)

ini <- mice(data, maxit = 0, printF = FALSE)

### Prediction matrix

# set the prediction matrix completely to 0, i.e. nothing predicts anything
pred <- ini$pred
pred[] <- 0

# set the variables that you want to impute as 1s in the predictor matrix so that they will be imputed 
pred[rownames(pred) %in% impvars, colnames(pred) %in% impvars] <- 1

# diagonal elements need to be 0s (i.e. one variable does not predict itself)
diag(pred) <- 0

### Imputation method

# get the method for imputation
meth <- ini$meth

# put to null "" for every combination that is no in the impvars 
meth[!names(meth) %in% impvars] <- "" 

### Order of imputation 
# the order in which variables are imputed matters in mice. 
# You can first put vars which can be predictive of the following vars 
visit <- ini$visit
visit <- visit[visit %in% impvars]
visit2 <- c("age_m_v2", "ethninfv3", "educm", "gsi", "msmoke", "can_m", "gestbir",
            "weight", "educm5", "income5", "f0300178", "sum_att_9m", "knicr_tbv_f9")

# run the imputation with m (number of datasets) and maxit (number of iterations)
## here 30 iterations and samples were set 
imp <- mice::mice(data, 
                  m = 30, 
                  maxit = 30, 
                  seed = 2020, 
                  pred = pred, 
                  meth = meth, 
                  visit = visit2)


saveRDS(imp, file.path(dir_data, "imp30-30_010321.rds"))
