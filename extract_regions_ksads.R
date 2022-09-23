library(data.table)
library(magrittr)
library(QDECR)

dir_wd <- "PATH/TO/YOUR/ABCD/WORKING/DIRECTORY"
source("postprocessing_functions.R")

FREESURFER_HOME <- Sys.getenv("FREESURFER_HOME")
SUBJECTS_DIR <- file.path(FREESURFER_HOME, "subjects")
FSAVERAGE <- file.path(SUBJECTS_DIR, "fsaverage")

#### OTHER VARS
hemis <- c("lh", "rh")
project_date <- "150722"
det <- "ksads_14_853_p"
out <- c("qdecr_area", "qdecr_volume", "qdecr_thickness")

#### INIT
setwd(dir_wd)

#### PATHS TO PROJECTS

out2 <- sub("qdecr_", "", out)
base_name <- paste0("_ksads_", project_date, ".", out2)

m <- 5 # number of models
paths_lh <- paths_rh <- list()

for (i in seq_along(out2)) {
  project_names <- expand.grid(".M", seq_len(m), base_name[i]) %>%
    apply(1, paste, collapse = "")
  
  projects_lh <- paste0("lh", project_names)
  projects_rh <- paste0("rh", project_names)
  
  paths_lh[[i]] <- file.path(dir_wd, "qdecr", projects_lh)
  paths_rh[[i]] <- file.path(dir_wd, "qdecr", projects_rh)
}

#### LOOP OVER ALL PROJECTS
# Let's loop over it and store the results:

lh_all <- unlist(paths_lh)
rh_all <- unlist(paths_rh)

lh_results <- lapply(lh_all, function(x) {
  vw <- qdecr_load(x)
  vw <- qdecr_update_path(vw, FREESURFER_HOME, SUBJECTS_DIR, dir_project = x)
  temp <- subset(summary(vw, annot = TRUE), variable == det)
  names(temp)[4] <- "mean_coef"
  
  if (nrow(temp) > 0) {
    model_number <- sub(".*lh.M+([0-9])+_.*", "\\1", x)
    model_outcome <- tools::file_ext(x)
    mm2 <- fread(paste(cmd = "grep -v '^#'", vw$stack$cluster.summary[[det]]))$V4
    cbind(outcome = model_outcome, model = model_number, hemi = "lh", mm2 = mm2, temp)
  } else {
    NA
  }
}) %>% do.call("rbind", .) %>%
  subset(!is.na(outcome))

rh_results <- lapply(rh_all, function(x) {
  vw <- qdecr_load(x)
  vw <- qdecr_update_path(vw, FREESURFER_HOME, SUBJECTS_DIR, dir_project = x)
  temp <- subset(summary(vw, annot = TRUE), variable == det)
  names(temp)[4] <- "mean_coef"
  
  if (nrow(temp) > 0) {
    model_number <- sub(".*rh.M+([0-9])+_.*", "\\1", x)
    model_outcome <- tools::file_ext(x)
    mm2 <- fread(paste(cmd = "grep -v '^#'", vw$stack$cluster.summary[[det]]))$V4
    cbind(outcome = model_outcome, model = model_number, hemi = "rh", mm2 = mm2, temp)
  } else {
    NA
  }
}) %>% do.call("rbind", .) %>%
  subset(!is.na(outcome))

all_results <- rbind(lh_results, rh_results)
write.table(all_results, "regions_extracted_from_ksads.txt", sep = "\t", quote = FALSE, row.names = FALSE)
