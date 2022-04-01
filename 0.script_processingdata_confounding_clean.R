###################
# 0. Data preparation
###################
# Project: brain morphology - ADHD, confounding
# Dataset: the Generation R Study @9
# Authors: Lorenza Dall'Aglio, as based on scripts by Sander Lamballais-Tessenshon
# NB previous steps: 
# data missings (missings (i.e. 888,999)) and missing system have been transformed into missing system in SPSS for readability in R


####
# STEP 1: Read in and merge the data
####

# load packages

x <- c("data.table", "foreign", "dplyr", "GGally", "car", "MASS", "mice")

lapply(x, require, character.only = T) 


# set directories - not specified here for privacy 

setwd("SETYOURDIRECTORY")
location <- "SETLOCATIONDATASETS"


# get the names of all datasets in the folder 

data <- list.files(pattern = ".sav")

data[1:3] <- data[3:1] # so the general data is first


# read in all files into a list

data2 <- list()

for (i in data){
  data2[[i]] <- read.spss(i, to.data.frame = T)
}



####
# STEP 2: Merge data
####

# some adjustments before merging

data3 <- lapply(data2, function(x) {names(x) <- tolower(names(x)); x}) # make all variables lowercase 

data4 <- lapply(data3, as.data.table) # all dataframes to data tables 


# Merge general data, smoking during pregnancy, drug use, BSI by idm (id of the mother) 
# (these datasets do not contain information on the IDC - id of the child, contrarily to other datasets)
# NB order before merging

# order

data4$`CHILD-ALLGENERALDATA_29012018.sav`[order(data4$`CHILD-ALLGENERALDATA_29012018.sav`$idm), ]

data4$`130815_ GEDRAGSGROEP Maternal smoking pregnancy.sav`[order(data4$`130815_ GEDRAGSGROEP Maternal smoking pregnancy.sav`$idm), ]

data4$`GR1003-BSI D1_22112016.sav`[order(data4$`GR1003-BSI D1_22112016.sav`$idm), ]


# merge

data4$temp <- merge(data4$`CHILD-ALLGENERALDATA_29012018.sav`, data4$`GR1003-BSI D1_22112016.sav`, by = intersect(names(data4$`CHILD-ALLGENERALDATA_29012018.sav`), names(data4$`GR1003-BSI D1_22112016.sav`)), all.x = T)

data4$general_maternal <- merge(data4$temp, data4$`130815_ GEDRAGSGROEP Maternal smoking pregnancy.sav`, by = intersect(names(data4$temp), names(data4$`130815_ GEDRAGSGROEP Maternal smoking pregnancy.sav`)), all.x = T)


# delete old datasets that have already been merged

data4$`CHILD-ALLGENERALDATA_29012018.sav` <- NULL

data4$`130815_ GEDRAGSGROEP Maternal smoking pregnancy.sav` <- NULL

data4$temp <- NULL



# Merge the remaining datasets by idc 

data4[1:7] <- data4[7:1] 

lapply(data4, function (x) {x[order(x$idc), ]}) # order all datasets by idc

dd <- Reduce(function(d1, d2) merge(d1, d2, by = intersect(names(d1), names(d2)), all.x = TRUE), data4) 
#reduce = member of apply family enabling to take into account result from previous iteration for the next

dd2 <- as.data.frame(dd)



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
          "age_m_v2", "educm5", "income5", "gestbir", "weight", "gestbir", 
          "msmoke", "can_m", "gsi", "f0300178", "educm","knicr_tbv_f9", "etiv_f9")


# filter the df for the variables we need (specified in vars)

dd3 <- dd2[ , names(dd2) %in% vars] 



####
# STEP 4: Variable type check 
####

str(dd3)

dd3$idc <- as.character(dd3$idc)
dd3$idm <- as.character(dd3$idm)
dd3$mother <- as.character(dd3$mother)



####
# STEP 4: Flow chart
#### 

# INCLUSION: Select individuals with inattention data & sufficient quality MRI data

# has ADHD data

dd4 <- dd3[!is.na(dd3$sum_att_9m), ] 


## has MRI consent 

dd5 <- dd4[dd4$mri_consent == "yes" & !is.na(dd4$mri_consent), ]


# EXCLUSION 

## scan type

dd6 <-  dd5[dd5$t1_asset_nii == "NO" & !is.na(dd5$t1_asset_nii), ]


## FS quality

dd6b <-  dd6[dd6$usable_fs_f9_final == "useable", ]


## Qdec

dd6c <- dd6b[dd6b$qdeclgi_f9 == "has data" & !is.na(dd6b$qdeclgi_f9), ]


## no braces 

dd7 <- dd6c[dd6c$braces_mri_f9 == "does not have braces", ]


## no incidental findings

dd8 <- dd7[dd7$exclude_incidental == "include" & !is.na(dd7$exclude_incidental), ]


## exclude twins

dd9 <- dd8[dd8$twin == "No", ]


## siblings (randomly exclude one of the siblings

set.seed(2020)

dd10 <- dd9[sample(nrow(dd9)), ] # randomly order data 

dd10 <- dd9[!duplicated(dd9$mother), ] # keep only one kid per mother


## Final sample = 2531 children 


####
# STEP 5: Collapse levels in categorical variables with many categories
####

#### facred: Factor level reduction; merges factor levels (mer) into a new level (new)
## Variables: 
#     x: A factor.
#     mer: Two or more levels that need to be merged.
#     new: The new level name. If NULL, the new name will be `paste(mer, collapse = clps)`.
#     clps: The `collapse` argument for `new`. This is only relevant if `new` is not null.



# function

facred <- function(x, mer, new = NULL, clps = " || "){
  le <- levels(x)
  le2 <- `names<-`(as.list(le), le)
  ind <- sapply(le2, function(y) y %in% mer) %>% unname
  if (is.null(new)) new <- paste(mer, collapse= clps) %>% as.character
  names(le2)[ind] <- new
  levels(x) <- le2
  x
}



# Check the initial levels of the variables

levels(dd10$ethninfv3)
levels(dd10$educm5)
levels(dd10$income5)
levels(dd10$educm)
levels(dd10$can_m)


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
    if (out %in% names(dd10)){
      old <- temp[[2]]
      new <- temp[[3]]
      dd10[,out] <- facred(dd10[,out], old, new)
    }
  }
}



####
# STEP 6: Descriptives & plots
####

# histograms 

png("hists.png")
par(mfrow = c(3,3))
hist(dd10$gsi, nclass = 100, main = "Maternal psychopathology (birth)")
hist(dd10$age_m_v2, nclass = 100, main = "Maternal age (birth)")
hist(dd10$gestbir, nclass = 100, main = "Gestational age (birth)")
hist(dd10$f0300178, nclass = 100, main = "Child IQ (age 5)")
hist(dd10$agechildbrainmrif9, nclass = 100, main = "Child age (at MRI)")
hist(dd10$sum_att_9m, nclass = 100, main = "Inattention (age 9)")
dev.off()


# Pairwise relations between the predictors

## cont by cont vars
png("scatterplots_covs.png")
par(mfrow=c(3,2))
plot(dd10$gsi, dd10$sum_att_9m, main="Maternal psychopathology and child inattention",
     xlab="Maternal psychopathology ", ylab="Child inattention")
plot(dd10$gsi, dd10$f0300178, main="Maternal psychopathology and child IQ",
     xlab="Maternal psychopathology", ylab="Child IQ")
plot(dd10$f0300178, dd10$sum_att_9m, main="Child IQ and inattention",
     xlab="Child IQ", ylab="Inattention")
plot(dd10$agechildbrainmrif9, dd10$sum_att_9m, main = "Child age and inattention", 
     xlab = "Child age", ylab = "Inattention")
dev.off()


# correlations

png("correlation_matrix.png")
ggcorr(dd10)
dev.off()

## cont by cat vars
png("boxplots_covs.png")
par(mfrow = c(3,2))
boxplot(gsi ~ educm5,data=dd10, xlab="Maternal education", ylab = "Maternal psychopathology", main = "Maternal psychopathology by education")
boxplot(sum_att_9m ~ adhdmed,data=dd10, xlab="Medication use", ylab = "Inattention", main = "Inattention symptoms by ADHD medication status")
boxplot(sum_att_9m ~ gender,data=dd10, xlab="Sex", ylab = "Inattention", main = "Inattention symptoms by sex")
boxplot(sum_att_9m ~ income5, data = dd10, xlab = "Income level", ylab = "Inattention", main = "Inattention symptoms by income")
boxplot(sum_att_9m ~ educm5, data = dd10, xlab = "Educational level", ylab = "Inattention", main = "Inattention symptoms by maternal educational level")
dev.off()


#### 
# STEP 8: OLR Assumption checks
####

# NB the total brain volume var is used as outcome as it is a proxy for surface area

# fit the model you want to test (with all confounders in) 
fit <- lm(dd10$knicr_tbv_f9 ~ dd10$sum_att_9m + dd10$gender + dd10$agechildbrainmrif9 + 
            dd10$ethninfv3 + dd10$educm5 + dd10$income5 + dd10$age_m_v2 + dd10$gsi +
            dd10$can_m + dd10$msmoke + dd10$f0300178)


# check linearity, normally-distributed residuals, homoscedasticity, outliers

png("diagnostic_plots_OLRassump.png")
par(mfrow = c(2,2))
plot(fit)
dev.off() 


# multicolinearity check

vif(fit) # variance inflation factors
sqrt(vif(fit)) > 2 


# autocorrelated errors test

durbinWatsonTest(fit)  


# global test of model assumptions

gvmodel <- gvlma(fit)

summary(gvmodel) 


# Normality of Residuals
# qq plot for studentized resid

## distribution of studentized residuals

sresid <- studres(fit)

#png("stud res.png")
hist(sresid, freq=FALSE,
     main="Distribution of Studentized Residuals", nclass = 100)
xfit<-seq(min(sresid),max(sresid),length=40)
yfit<-dnorm(xfit)
lines(xfit, yfit)
#dev.off()



## Evaluate Nonlinearity
# component + residual plot

crPlots(fit)

# Ceres plots

ceresPlots(fit)



####
# STEP 7: Prepare data for imputation
####

# Get the N & % missingness

missings_n <-  apply(dd10, 2, function(col) {sum(is.na(col))})

write.csv(missings_n,"PATH/missings_N.csv", row.names = F)

missings_percent <- apply(dd10, 2, function(col) {((sum(is.na(col)))/length(col)*100)})

write.csv(missings_percent,"PATH/missing_stats_percent.csv", row.names = F)


# Get summary stats

summarystats <- summary(dd10) 

write.csv(summarystats,"PATH/descriptive_stats.csv", row.names = F)


# Get the sd for Suppl T 2

gsisd <- sd(dd10$gsi, na.rm = T)

age_m_v2_sd <- sd(dd10$age_m_v2, na.rm = T)

gestbir_sd <- sd(dd10$gestbir, na.rm = T) 

iq_sd <- sd(dd10$f0300178, na.rm = T)

age_sd <- sd(dd10$agechildbrainmrif9, na.rm = T)


sd_all <- cbind(gsisd, age_m_v2_sd, gestbir_sd, iq_sd, age_sd, agg_sd)  

write.csv(sd_all,"PATH/sds.csv", row.names = F)



####
# STEP 8: impute
####


#dryrun mice (page 35 mice guide)

ini <- mice(dd10, maxit = 0, printF = FALSE)
ini


### Variables 

# specify the variables that can aid the prediction (i.e. help predict the missingness in other variables)

predictors_for_imputation <- c("age_m_v2", "educm", "income5", "msmoke", "can_m", "sum_att_9m", "knicr_tbv_f9")


# specify the variables that you want to impute. It's important to add in the predictors for imputation 
## this is important because if our predictors have NAs, our output of the imputation will also present NAs 

impvars <- c("gestbir", "weight", "educm5", "f0300178", "gsi", "ethninfv3", "income5", "msmoke", predictors_for_imputation)  


### Prediction matrix

# set the prediction matrix completely to 0, i.e. nothing predicts anything

pred <- ini$pred
pred[] <- 0


# set the variables that you want to impute as 1s in the predictor matrix so that they will be imputed 

pred[rownames(pred) %in% impvars, colnames(pred) %in% impvars] <- 1

# diagonal elements need to be 0s (i.e. one variable does not predict itself)

diag(pred) <- 0

pred #check prediction matrix


### Imputation method

# get the method for imputation

meth <- ini$meth


# put to null "" for every combination that is no in the impvars 

meth[!names(meth) %in% impvars] <- "" 

meth


### Order of imputation 
# the order in which variables are imputed matters in mice. 
# You can first put vars which can be predictive of the following vars 

visit <- ini$visit

visit <- visit[visit %in% impvars]

visit2 <- c("age_m_v2", "ethninfv3", "educm", "gsi", "msmoke", "can_m", "gestbir",
            "weight", "educm5", "income5", "f0300178", "sum_att_9m", "knicr_tbv_f9")



# run the imputation with m (number of datasets) and maxit (number of iterations)
## here 30 iterations and samples were set 

dd11 <- mice::mice(dd10, 
                  m = 30, 
                  maxit = 30, 
                  seed = 2020, 
                  pred = pred, 
                  meth = meth, 
                  visit = visit2)
				  

saveRDS(dd11, "PATH/imp30-30_010321.rds")

dd11 <- readRDS("imp30-30_010321.rds") 



####
# Step 9: Imputation diagnostics
####

# Checking for convergence (pg. 39 mice guide) 
## convergence occurs when lines overlap, but no general trend is present 

png("imputation_diagnostics_plot1.png")
plot(dd11, c("gestbir", "weight", "educm5"))
dev.off()

png("imputation_diagnostics_plot2.png")
plot(dd11, c("f0300178", "gsi", "ethninfv3"))
dev.off()

png("imputation_diagnostics_plot3.png")
plot(dd11, c("income5", "msmoke", "can_m"))
dev.off()



# Checking for whether the imputed values resemble the observed values 
## for large samples density plots are used (for continuous vars)  

densityplot(dd11)



