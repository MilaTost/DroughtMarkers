rm(list=ls())
# BeechProj: Calculate Ghat ----------------------------------------------------
### Load packages --------------------------------------------------------------
library(stringr)
library(data.table)
library(doMC)
library(parallel)
library(foreach)
library(rrBLUP)
library(R.utils)
### Set the parameters ---------------------------------------------------------
cores <- 4
time_0 <- Sys.time()
VCF_file_dir <- "Beech_project/Romanian_Beech_Data/Ghat_test_trials/Res/"
working_dir <- "Beech_project/Romanian_Beech_Data/Ghat_test_trials/"
result_dir <- "Beech_project/Romanian_Beech_Data/Ghat_test_trials/Res/"
input_file_dir <- "Beech_project/Romanian_Beech_Data/Phenotypic_data/"
VCF_file_dir <- "Beech_project/Romanian_Beech_Data/Genotypic_data/"
VCF_file <- "all_merged.norm.filtered.vcf"
VCF_file_describ <- "RomanianBeech"
### Load the data sets ---------------------------------------------------------
#### Load the VCF file
cat("In this script NAs are replaced by the mean.", "\n", "\n")
run_description <- "100Kperm_noSolomon_outliers_NOT_removed_MAFfilt_miss_filt"
cat("The data set loading starts","\n")
VCF <- fread(paste0(VCF_file_dir, VCF_file), skip = 110)
VCF <- as.data.frame(VCF)
SNP_IDs <- paste0(VCF[,1], "_", VCF[,2])
colnames(VCF)[1] <- "CHROM"
time_1 <- Sys.time()
cat("The data set was loaded after:","\n")
print(time_1 - time_0)
#### Load the other input files
##### Load phenotypic and environmental data
setwd(input_file_dir)
pheno <- read.table("2024_06_24_RomanianBeech_all_traits.txt", sep = "\t", header = TRUE)
locations <- read.table("2023_09_01_RoBe_Location_and_weather_data.txt", sep = "\t", header = TRUE)
env_PC <- read.table("2025_Env_Covar.txt", sep = "\t", header = TRUE)
setwd(working_dir)
### Get gt file ----------------------------------------------------------------
### Change colnames of VCF file ------------------------------------------------
sample_names <- colnames(VCF)[10:length(colnames(VCF))]
# something is wrong with the sample names
split_sample_names <- str_split(sample_names, "-")
split_sample_names <- matrix(data = unlist(split_sample_names),
                             nrow = length(colnames(VCF)[10:length(colnames(VCF))]),
                             ncol = 4, byrow = TRUE)
new_sample_names <- split_sample_names[, 2]
index_Ruia <- which(str_sub(new_sample_names, 1, 3) == "RuB")
index_Tampa <- which(str_sub(new_sample_names, 1, 3) == "TaB")
index_Solomon <- which(str_sub(new_sample_names, 1, 3) == "SoB")
index_Sinpetru <- which(str_sub(new_sample_names, 1, 3) == "SiB")
index_Lupului <- which(str_sub(new_sample_names, 1, 3) == "LuB")
new_sample_names[index_Ruia] <- str_replace(new_sample_names[index_Ruia], "RuB", "Ruia")
new_sample_names[index_Tampa] <- str_replace(new_sample_names[index_Tampa], "TaB", "Tampa")
new_sample_names[index_Solomon] <- str_replace(new_sample_names[index_Solomon], "SoB", "Solomon")
new_sample_names[index_Sinpetru] <- str_replace(new_sample_names[index_Sinpetru], "SiB", "Sinpetru")
new_sample_names[index_Lupului] <- str_replace(new_sample_names[index_Lupului], "LuB", "Lupului")
colnames(VCF)[10:length(colnames(VCF))] <- new_sample_names
cat("How many Solomon 53 samples are in the data set:", length(which(colnames(VCF) == "Solomon53")), "\n")
cat("The column names of the VCF file were changed!", "\n")
#new_VCF <- VCF[,-which(colnames(VCF) == "Solomon53")[2]]
#VCF <- new_VCF
new_VCF <- VCF[,which(colnames(VCF) != "Solomon34")]
cat("How many Solomon 53 samples are in the data set:", length(which(colnames(new_VCF) == "Solomon53")), "\n")
VCF <- new_VCF
### Extract genotypes function -------------------------------------------------
cat("VCF file contains now:","\n",
    nrow(VCF),"markers and", length(colnames(VCF)[10:length(colnames(VCF))]), "individuals","\n","\n")
cat("Genotypes are getting extracted.", "\n", "\n")
extract_genotypes <- function(VCF){
  gt_samples <- VCF[, 10:ncol(VCF)]
  gt_samples <- as.matrix(gt_samples)
  nr_of_markers_in_subset <- nrow(gt_samples)
  cat("Genotypes are getting extracted.", "\n")
  cat("Extract genotypes per SNP!", "\n")
  run_extract_gt <- function(i){
    new_geno <- vector(length = ncol(gt_samples))
    new_geno[1:length(new_geno)] <- 9
    index_0 <- which(str_detect(gt_samples[i,], "^0/0"))
    index_2 <- which(str_detect(gt_samples[i,], "^1/1"))
    index_het1 <- which(str_detect(gt_samples[i,], "^1/0"))
    index_het2 <- which(str_detect(gt_samples[i,], "^0/1"))
    new_geno[index_0] <- 0
    new_geno[index_2] <- 2
    new_geno[index_het1] <- 1
    new_geno[index_het2] <- 1
    return(new_geno)
  }
  gt <- mclapply(1:nr_of_markers_in_subset, run_extract_gt)
  gt <- unlist(gt)
  gt <- matrix(data = gt, nrow = nr_of_markers_in_subset,
               ncol = ncol(gt_samples), byrow = TRUE)
  cat("Extraction of", nrow(gt),"genotypes is finished.", "\n", "\n")
  return(gt)
}
gt <- extract_genotypes(VCF = VCF)
gt <- as.data.frame(gt)
colnames(gt) <- colnames(VCF)[10:ncol(VCF)]
rownames(gt) <- SNP_IDs
## pheno to phenotypes
### Phenotypic data ------------------------------------------------------------
traits <- colnames(pheno)[4:length(colnames(pheno))]
phenotypes <- t(pheno[, 4:length(colnames(pheno))])
colnames(phenotypes) <- pheno$Sample_ID
phenotypes <- as.data.frame(phenotypes)
phenotypes <- cbind(phenotypes[3:13, ])
### Filtering of data sets -----------------------------------------------------
#### Filter location data
filter_for_available_samples_in_location_data <- function(i){

  if(any(colnames(phenotypes)[i] == env_PC$FID)){
    index <- which(colnames(phenotypes)[i] == env_PC$FID)
    return(index)}
}
index_loc_data <- mclapply(1:length(phenotypes), filter_for_available_samples_in_location_data)
index_loc_data <- unlist(index_loc_data)
new_env_PC <- env_PC[index_loc_data, ]
env_PC <- new_env_PC
cat("The phenotype file contains", length(colnames(phenotypes)), "individuals.", "\n",
    "The genotype file contains", length(colnames(gt)), "individuals.", "\n",
    "The location file contains", nrow(env_PC), "individuals.", "\n")
cat("The filtering of phenotypic and genotypic file was done!", "\n", "\n")
#### Filter genotypes
filter_for_available_samples_in_gt <- function(i){
  if(any(env_PC$FID[i] == colnames(gt))){
    index <- which(env_PC$FID[i] == colnames(gt))
    return(index)}
}
index_gen <- mclapply(1:nrow(env_PC), filter_for_available_samples_in_gt)
index_gen <- unlist(index_gen)
new_gen <- gt[, index_gen]
rownames(new_gen) <- SNP_IDs
gt <- new_gen
rm(new_gen)
cat("The phenotype file contains", length(colnames(phenotypes)), "individuals.", "\n",
    "The genotype file contains", length(colnames(gt)), "individuals.", "\n")
#### Filter phenotypes
colnames_gt <- colnames(gt)
filter_for_available_samples_in_phen <- function(i){
  if(any(colnames_gt[i] == colnames(phenotypes))){
    index <- which(colnames_gt[i] == colnames(phenotypes))
    return(index)}
}
index_phen <- mclapply(1:length(colnames_gt), filter_for_available_samples_in_phen)
index_phen <- unlist(index_phen)
new_phen <- phenotypes[, index_phen]
phenotypes <- new_phen
#### Filter location data
filter_for_available_samples_in_location_data <- function(i){
  if(any(colnames(phenotypes)[i] == env_PC$FID)){
    index <- which(colnames(phenotypes)[i] == env_PC$FID)
    return(index)}
}
index_loc_data <- mclapply(1:length(phenotypes), filter_for_available_samples_in_location_data)
index_loc_data <- unlist(index_loc_data)
new_env_PC <- env_PC[index_loc_data, ]
env_PC <- new_env_PC
cat("The phenotype file contains", length(colnames(phenotypes)), "individuals.", "\n",
    "The genotype file contains", length(colnames(gt)), "individuals.", "\n",
    "The location file contains", nrow(env_PC), "individuals.", "\n")
### Remove markers with high missingness ---------------------------------------
cat("Calculate missingness at markers!", "\n")
calc_missingness <- function(dt){
  calc_missingness_per_marker <- function(i){
    nr_missing_obs <- length(which(dt[i ,] == 9))
    rate_of_missing_obs <- nr_missing_obs/nrow(dt)
  }
  missingness_list <- mclapply(1:nrow(dt), calc_missingness_per_marker)
  missingness_dt <- unlist(missingness_list)
  missingness_dt <- matrix(data = missingness_dt, ncol = 1, nrow = nrow(dt))
  missingness_dt <- as.data.frame(missingness_dt)
  rownames(missingness_dt) <- rownames(dt)
  colnames(missingness_dt) <- "Rate_of_missingness"
  return(missingness_dt)
}
dt_missingness <- calc_missingness(dt = gt)
write.table(dt_missingness, paste0(result_dir, "Missingness_per_marker_before_filtering_", run_description, ".txt"))
cat("Remove markers with too high missing markers!", "\n")
missingness_value <- 0.05
remove_markers_to_high_missingness <- function(Y){
  replace_per_marker <- function(i){
    if(all(Y[i ,] == 9)){
      index <- NA
      return(index)
    }
    if(length(which(Y[i ,] == 9)) >= round((ncol(Y)*missingness_value), 0)){
      index <- NA
      return(index)
    }
    if(length(which(Y[i ,] == 9)) <= round((ncol(Y)*missingness_value), 0)){
      index <- i
      return(index)
    }
  }
  index_mat <- mclapply(1:nrow(Y), replace_per_marker)
  index_mat <- unlist(index_mat)
  index <- which(!is.na(index_mat))
  Y <- Y[index, ]
  return(Y)
}
new_gt <- remove_markers_to_high_missingness(Y = gt)
cat("Number of missing observations in the data set is:", length(which(new_gt==9)),"\n")
cat("Before we had", nrow(gt),"\n")
cat("Now we have", nrow(new_gt),"\n")
gt <- new_gt
rm(new_gt)
cat("Replace missing observations!", "\n")
SNP_IDs <- rownames(gt)
### Remove markers with high missingness ---------------------------------------
remove_markers_too_high_MAF <- function(Y){
  replace_per_marker <- function(i){
    tot_obs <- 2*length(which(Y[i,]!=9))
    af_freq_HOM_Ref <- length(which(Y[i,]==0))
    af_freq_HET <- length(which(Y[i,]==1))
    af_freq_HOM_alt <- length(which(Y[i,]==2))
    af_MAF <- (af_freq_HET+af_freq_HOM_alt)/tot_obs
    if(af_MAF <= 0.01){
      index <- NA
    }
    if(af_MAF > 0.01){
      index <- i
    }
    return(index)
  }
  index_mat <- mclapply(1:nrow(Y), replace_per_marker)
  index_mat <- unlist(index_mat)
  index <- which(!is.na(index_mat))
  SNP_IDs <- rownames(Y)[index]
  Y <- Y[index, ]
  rownames(Y) <- SNP_IDs
  return(Y)
}
new_gt <- remove_markers_too_high_MAF(Y = gt)
cat("Number of markers before filtering for MAF:", nrow(gt), "\n",
    "Number of markers after filtering for MAF:", nrow(new_gt), "\n")
calc_MAF <- function(Y){
  replace_per_marker <- function(i){
    tot_obs <- 2*length(which(Y[i,]!=9))
    af_freq_HOM_Ref <- length(which(Y[i,]==0))
    af_freq_HET <- length(which(Y[i,]==1))
    af_freq_HOM_alt <- length(which(Y[i,]==2))
    af_MAF <- (af_freq_HET+af_freq_HOM_alt)/tot_obs
    return(af_MAF)
  }
  dt_MAF <- foreach(i = 1:nrow(Y), .combine=rbind) %do% replace_per_marker(i)
  dt_MAF <- as.data.frame(dt_MAF)
  return(dt_MAF)
}
dt_MAF_new_gt <- calc_MAF(Y = new_gt)
gt <- new_gt
SNP_IDs <- rownames(new_gt)
### Add randomly generated traits ----------------------------------------------
# how many traits are necessary? 100 traits?
### Group division -------------------------------------------------------------
pop1 <- "Ruia"
pop3 <- "Lupului"
pop2 <- "Solomon"
pop4 <- "Tampa"
pop5 <- "Sinpetru"
group1 <- c(pop1, pop3)
group2 <- c(pop4, pop5)
index_pop1 <- which(str_detect(colnames(gt), pop1))
index_pop3 <- which(str_detect(colnames(gt), pop3))
index_pop4 <- which(str_detect(colnames(gt), pop4))
index_pop5 <- which(str_detect(colnames(gt), pop5))
gt_group_1 <- gt[ ,c(index_pop1, index_pop3)]
gt_group_2 <- gt[ ,c(index_pop4, index_pop5)]
index_pop1 <- which(str_detect(phenotypes$Pop, pop1))
index_pop3 <- which(str_detect(phenotypes$Pop, pop3))
index_pop4 <- which(str_detect(phenotypes$Pop, pop4))
index_pop5 <- which(str_detect(phenotypes$Pop, pop5))
pheno_group_1 <- phenotypes[, c(index_pop1, index_pop3)]
pheno_group_2 <- phenotypes[, c(index_pop4, index_pop5)]
cat("VCF file contains now:","\n",
    nrow(gt),"markers and", length(colnames(gt)), "individuals","\n","\n")
cat("Genotypes are getting extracted.", "\n", "\n")
### Replace missing values -----------------------------------------------------
# do the replacement first based on the mean and use later the imputed data set
cat("Replacing of missing values by population mean starts!", "\n")
# in case there are still missing observations, we replace them by the marker mean
replace_missing_markers_by_mean <- function(i){
  gt[,which(gt[i, ]==9)]<- NA
  marker_gt <- gt[i, ]
  if(any(is.na(marker_gt))){
    cat("Missing observations found.", "\n")
    marker_gt[which(is.na(marker_gt))] <- mean(as.numeric(marker_gt), na.rm = TRUE)
    marker_gt <- as.numeric(marker_gt)
    return(marker_gt)
  }
  if(all(!is.na(marker_gt))){
    #cat("No missing observations.", "\n")
    marker_gt <- as.numeric(marker_gt)
    return(marker_gt)
  }
}
new_gt <- mclapply(1:nrow(gt), replace_missing_markers_by_mean)
cat("Replacing of missing values was done!", "\n")
new_gt <- matrix(data = unlist(new_gt),
                 nrow = nrow(gt),
                 ncol = ncol(gt),
                 byrow = TRUE)
colnames(new_gt) <- colnames(gt)
rownames(new_gt) <- SNP_IDs
new_gt <- as.data.frame(new_gt)
gt <- new_gt
#cat("The data set looks like this:", str(gt)[1:100], "\n")
rm(new_gt)
gt <- as.matrix(gt)
### Phenotypic data ------------------------------------------------------------
traits <- rownames(phenotypes)
Sample_names <- colnames(phenotypes)
phenotypes <- as.data.frame(phenotypes)
new_env_PC <- t(env_PC[3:ncol(env_PC)])
colnames(new_env_PC) <- env_PC[,2]
number_of_effective_markers <- 7000
time_1 <- Sys.time()
### Load the parameters for the Ghat function ----------------------------------
cat("The parameters are getting loaded.", "\n")
NA_replacement <- "by population mean"
calc_change <- "between groups"
grouping <- "based on measured gsp"
other_variables_in_model <- "subpop"
group <- "Mean_percip"
group1 <- c("Ruia", "Solomon")
group2 <- c("Tampa", "Sinpetru")
permutations_to_run <- 100000
cat("This script is using the following parameters:","\n",
    "NAs are replaced:","\t","\t", NA_replacement, "\n",
    "Change is calculated as:","\t", calc_change, "\n",
    "Other variables in the model:", "\t", other_variables_in_model, "\n",
    "Grouping is done based on:", "\t", grouping, "\n",
    "\t","\t","\t","\t", group, "\n",
    "Number of permutations:", "\t", permutations_to_run,"\n",
    "Number_of_effective_markers:","\t", number_of_effective_markers,"\n",
    "\n")
### Group identification -------------------------------------------------------
### Randomly sample the genotypes file -----------------------------------------
#### Based on the genotypic file
cat("Group 1 consists out of", ncol(gt_group_1), "different populations which are", group1, "\n",
    "Group 2 consists out of", ncol(gt_group_2), "different populations which are", group2, "\n")
### Get only one phenotype -----------------------------------------------------
pheno <- phenotypes
traits <- rownames(pheno)
get_only_one_phenotype <- function(i){
  pheno_name <- rownames(pheno)[i]
  pheno_dt <- t(pheno[i,])
  colnames(pheno_dt) <- pheno_name
  return(pheno_dt)
}
dt_phenotypes <- mclapply(1:nrow(pheno), get_only_one_phenotype)
n_phenotypes <- length(lengths(dt_phenotypes))
cat("In total:", n_phenotypes, "traits are tested!", "\n")
### Calculate allele frequency change between groups ---------------------------
calculate_change_in_allele_freq <- function(group_1,
                                            group_2){
  tot_alleles_group2 <- rowSums(group_2, na.rm = TRUE)
  tot_alleles_group1 <- rowSums(group_1, na.rm = TRUE)
  nr_alleles_group2 <- 2*ncol(group_2)
  nr_alleles_group1 <- 2*ncol(group_1)
  change <- (tot_alleles_group2/nr_alleles_group2) - (tot_alleles_group1/nr_alleles_group1)
  return(change)
}
change <- calculate_change_in_allele_freq(group_1 = gt_group_1,
                                          group_2 = gt_group_2)
cat("The real allele frequency difference was calculated.","\n")
change <- as.numeric(change)
cat("The real allele frequency difference was calculated.","\n")
change <- as.numeric(change)
dt_change <- cbind(SNP_IDs, change)
dt_change <- as.data.frame(dt_change)
colnames(dt_change) <- c("SNP_ID", "Allele_freq_difference")
write.table(dt_change, paste0(result_dir, Sys.Date(),"_Allele_freq_difference_", run_description, ".txt"),
            sep = "\t",
            row.names = FALSE)
# only use meaningful predictions, calculate p-values
# add this part in another script
cat("The allele frequency change as a gradient slope was calculated.","\n")
### Estimate marker effects ----------------------------------------------------
cat("Test the estimation of marker effects.","\n")
genotypes <- gt
phenotypes <- dt_phenotypes
location_data <- locations
scaled_grad <- group
populations <- gsub("[^a-zA-Z]", "", colnames(genotypes))
gt_populations <- levels(as.factor(populations))
### Calculate marker effects for all traist ------------------------------------
cat("The estimation of marker effects starts.","\n")
estimate_marker_effects <- function(genotypes,
                                    phenotypes,
                                    location_data,
                                    scaled_grad){
  populations <- gsub("[^a-zA-Z]", "", colnames(genotypes))
  gt_populations <- levels(as.factor(populations))
  new_location_data <- aggregate(location_data[, c(4,6:18)], list(location_data[, 1]), mean)
  estimate_marker_effects_per_pheno <- function(i){
    cat("The marker effects are getting calculated for",
        colnames(phenotypes[[i]]),"\n")
    rownames(genotypes) <- rownames(genotypes)
    new_genotypes <- genotypes[,1:ncol(genotypes)]
    genotypes <- new_genotypes
    pheno_data <- phenotypes[[i]]
    pheno_unlisted <- unlist(pheno_data)
    pheno_data <- as.data.frame(pheno_unlisted)
    colnames(pheno_data) <- colnames(phenotypes[[i]])
    new_pheno_data <- as.data.frame(cbind(rownames(pheno_data), pheno_data))
    #pheno_data <- new_pheno_data[order(new_pheno_data[, 1]), ]
    pheno_data <- new_pheno_data
    colnames(pheno_data) <- c("Sample", colnames(phenotypes[[i]]))
    pheno_data$populations <- gsub("[^a-zA-Z]", "", rownames(pheno_data))
    env_PC$Pop <- gsub("[^a-zA-Z]", "", env_PC$FID)
    # For 2020
    if(str_detect(colnames(phenotypes[[i]]), "2020")){
      create_info_tab <- function(a){
      new_info <- env_PC[which(gt_populations[a] == env_PC$Pop), c(2:4)]
      new_envPC <- env_PC[which(gt_populations[a] == env_PC$Pop), c(5)]
      new_dt <- cbind(new_info, new_envPC)
      return(new_dt)
      }
    }
    # For 2021
    if(str_detect(colnames(phenotypes[[i]]), "2021")){
    create_info_tab <- function(a){
      new_info <- env_PC[which(gt_populations[a] == env_PC$Pop), c(2:4)]
      new_envPC <- env_PC[which(gt_populations[a] == env_PC$Pop), c(6)]
      new_dt <- cbind(new_info, new_envPC)
      return(new_dt)
    }
      }
    # For 2022
    if(str_detect(colnames(phenotypes[[i]]), "2022")){
    create_info_tab <- function(a){
      new_info <- env_PC[which(gt_populations[a] == env_PC$Pop), c(2:4)]
      new_envPC <- env_PC[which(gt_populations[a] == env_PC$Pop), c(7)]
      new_dt <- cbind(new_info, new_envPC)
      return(new_dt)
      }
    }
    new_info_dt <- foreach(a = 1:length(populations), .combine = rbind) %do% create_info_tab(a)
    pheno_dt <- cbind(pheno_data[,1:3], new_info_dt)
    cat("The column names from the group and pheno dt are matching:", all(match(gsub("[^a-zA-Z]", "", pheno_dt[,1]), pheno_dt[,3])), "\n")
    pheno_dt <- as.data.frame(pheno_dt)
    colnames(pheno_dt) <- c("Sample","Trait","Population","ID", "spet_PC1", "spet_PC2", "EnvPC1_of_year")
    pheno_dt$Population <- as.factor(pheno_dt$Population)
    pheno_dt$EnvPC1_of_year <- as.factor(pheno_dt$EnvPC1_of_year)
    pheno_dt$spet_PC1 <- as.factor(pheno_dt$spet_PC1)
    pheno_dt$spet_PC2 <- as.factor(pheno_dt$spet_PC2)
    # always check if the names of the genotypes are included or not
    # transpose the genotype matrix
    gt <- t(genotypes)
    # Remove the samples with missing observations from the data sets
    # Exclude traits with only 0s and NAs
    pheno_dt <- pheno_dt[which(is.na(pheno_dt[,2])==FALSE),]
    gt <- gt[which(is.na(pheno_dt[,2])==FALSE),]
    if(all(pheno_dt[,2] == 0)){
      effects <- NA
    }
    if(any(pheno_dt[,2] != 0)){
      result <- mixed.solve(pheno_dt[,2],
                            Z= as.matrix(gt),
                            X=model.matrix(pheno_dt[,2]~pheno_dt[,7]),
                            SE=FALSE, return.Hinv=FALSE,
                            method="ML")
      effects <- result$u
    }
    return(effects)
  }
  marker_effect_estimates <- mclapply(1:length(dt_phenotypes), estimate_marker_effects_per_pheno)
  return(marker_effect_estimates)
}
marker_effect_estimates <- estimate_marker_effects(genotypes = gt,
                                                   phenotypes = dt_phenotypes,
                                                   location_data = locations,
                                                   scaled_grad = group)
dt_SNP_effects_D132020 <- matrix(data = unlist(marker_effect_estimates[[2]]),
                                nrow = length(SNP_IDs), ncol = 1)
dt_SNP_effects_D132021 <- matrix(data = unlist(marker_effect_estimates[[6]]),
                                nrow = length(SNP_IDs), ncol = 1)
dt_SNP_effects_D132022 <- matrix(data = unlist(marker_effect_estimates[[10]]),
                                nrow = length(SNP_IDs), ncol = 1)
dt_SNP_effects <- rbind(cbind(rep("Delta13_2020", length(SNP_IDs)), SNP_IDs, dt_SNP_effects_D132020),
                        cbind(rep("Delta13_2021", length(SNP_IDs)), SNP_IDs, dt_SNP_effects_D132021),
                        cbind(rep("Delta13_2022", length(SNP_IDs)), SNP_IDs, dt_SNP_effects_D132022))
dt_SNP_effects <- as.data.frame(dt_SNP_effects)
colnames(dt_SNP_effects) <- c("Trait", "SNP_ID", "Marker_effect_estimates")
write.table(dt_SNP_effects, paste0(result_dir, Sys.Date(),"_Marker_effect_estimates_", run_description, ".txt"),
            sep = "\t",
            row.names = FALSE)
#marker_effect_estimates <- marker_effect_estimates[which(!is.na(marker_effect_estimates))]
#n_phenotypes <- length(lengths(marker_effect_estimates))
#is.na(marker_effect_estimates[[6]])
cat("The marker effect estimates were calculated.","\n")
### Calculate Ghat -------------------------------------------------------------
#source("C:/Users/mtost/Documents/Ghat_maize_diversity_panel_project/New_versions_of_Ghat_function/2023_10_24_new_Ghat.R")
source("/usr/users/mtost/Ghat_Maize_diversity_project/New_versions_of_Ghat/2022_07_11_Ghat_new_version.R")
cat("Ghat calculation starts.","\n")
run_Ghat_test <- function(i){
  #cat("Ghat is getting calculated for",
  #names(marker_effect_estimates)[[i]],"\n")
  if(all(!is.na(marker_effect_estimates[[i]]))){
    new_marker_effect_estimates <- marker_effect_estimates[[i]]
    Ghat_results <- new_Ghat(effects = new_marker_effect_estimates,
                             change = change,
                             method = "scale",
                             perms = permutations_to_run,
                             plot = F,
                             num_eff = number_of_effective_markers)
    result <- c(Ghat_results$Ghat,Ghat_results$Cor,Ghat_results$p.val)
    #cat(result,"\n")
  }
  if(all(is.na(marker_effect_estimates[[i]]))){
    result <- c(NA, NA, NA)
    #cat(result,"\n")
  }
  return(result)
}
#length(marker_effect_estimates)
number_of_effective_markers <- number_of_effective_markers
dt_Ghat_per_phenotype_real <- mclapply(1:length(marker_effect_estimates), run_Ghat_test)
dt_Ghat_real <- matrix(data = unlist(dt_Ghat_per_phenotype_real),
                       nrow = as.numeric(length(unlist(dt_Ghat_per_phenotype_real))/3),
                       ncol = 3, byrow = TRUE)
dt_Ghat_real <- as.data.frame(dt_Ghat_real)
colnames(dt_Ghat_real) <- c("Ghat","Correlation_coefficient","P_value")
dt_Ghat_real$Correlation_coefficient <- round(dt_Ghat_real$Correlation_coefficient, 4)
dt_Ghat_real$P_value <- round(dt_Ghat_real$P_value, 4)
rownames(dt_Ghat_real) <- traits[1:(length(unlist(dt_Ghat_per_phenotype_real))/3)]
cat("When we used", number_of_effective_markers, "markers:", "\n")
p_values_real_traits <- dt_Ghat_real$P_value
cat(length(which(p_values_real_traits < 0.05)), "real traits were detected to be significant", "\n")
cat("The result data set is getting created now!","\n")
dt_Ghat_real <- cbind(traits, dt_Ghat_real)
colnames(dt_Ghat_real) <- c("Trait","Ghat","Correlation_coefficient","P_value")
cat("The result data set is getting created now!","\n")
write.table(dt_Ghat_real, paste0(result_dir, Sys.Date(),"_Ghat_outliers_NOT_removed_all_traits_", run_description,".txt"),
            sep = "\t",
            row.names = FALSE)
cat("Wohooooo!!!!", "\n")
q()
