rm(list=ls())
# BeechProj: EAA ---------------------------------------------------------------
### Load packages --------------------------------------------------------------
library(stringr)
library(data.table)
library(parallel)
library(foreach)
library(rrBLUP)
library(R.utils)
library(doParallel)
# install.packages("devtools")
#devtools::install_github("bcm-uga/lfmm")
#devtools::install_github("bcm-uga/LEA")
library(lfmm)
library(LEA)
### Set the parameters ---------------------------------------------------------
cores <- 4
args <- commandArgs(trailingOnly=T)
i <- args[1]
i <- as.numeric(i)
latent_factor_K <- 3
latent_factor_K <- as.numeric(latent_factor_K)
lambda_value <- 1e-10
method <- "LFMM2"
run_description <- paste0("K", "_", latent_factor_K, "_", method)
### Load files -----------------------------------------------------------------
VCF_file <- "all_merged.norm.filtered.vcf"
VCF_file_describ <- "RomanianBeech"
VCF_file_dir <- "/Beech_project/Romanian_Beech_Data/Genotypic_data/"
working_dir <- "/Beech_project/Romanian_Beech_Data/EAA/LFMM/LFMM2_final_run/"
result_dir <- "/Beech_project/Romanian_Beech_Data/EAA/LFMM/LFMM2_final_run/Results/"
input_file_dir <- "/Beech_project/Romanian_Beech_Data/Phenotypic_data/"
lfmm_input_dir <- "/Beech_project/Romanian_Beech_Data/EAA/LFMM/LFMM2_final_run/new_input_dir/"
plot_dir <- "/Beech_project/Romanian_Beech_Data/EAA/LFMM/LFMM2_final_run/Plot_dir/"
### Load the data sets ---------------------------------------------------------
#### Load the VCF file
cat("In this script NAs are replaced by the impute function of the LEA package.", "\n", "\n")
cat("The data set loading starts","\n")
VCF <- fread(paste0(VCF_file_dir, VCF_file), skip = 110)
VCF <- as.data.frame(VCF)
SNP_IDs <- paste0(VCF[,1], "_",VCF[,2])
#colnames(VCF)[1] <- "CHROM"
#### Load the other input files
##### Load phenotypic and environmental data
setwd(input_file_dir)
pheno <- read.table("2024_06_24_RomanianBeech_all_traits.txt", sep = "\t", header = TRUE)
locations <- read.table("2024_06_19_All_environmental_variables.txt", sep = "\t", header = TRUE)
env_PCs <- read.table("env_PCs.txt", sep = "\t", header = TRUE)
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
new_VCF <- VCF[,-which(colnames(VCF) == "Solomon53")[2]]
#new_VCF
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
rownames(gt) <- rownames(VCF)
### Phenotypic data ------------------------------------------------------------
traits <- colnames(pheno)[4:length(colnames(pheno))]
phenotypes <- t(pheno[, 4:length(colnames(pheno))])
colnames(phenotypes) <- pheno$Sample_ID
phenotypes <- as.data.frame(phenotypes)
### Filtering of data sets -----------------------------------------------------
#### Filter location data
filter_for_available_samples_in_location_data <- function(i){

  if(any(colnames(phenotypes)[i] == locations$Sample_ID)){
    index <- which(colnames(phenotypes)[i] == locations$Sample_ID)
    return(index)}
}
index_loc_data <- mclapply(1:length(phenotypes), filter_for_available_samples_in_location_data, mc.cores = cores)
index_loc_data <- unlist(index_loc_data)
new_loc <- locations[index_loc_data, ]
locations <- new_loc
cat("The phenotype file contains", length(colnames(phenotypes)), "individuals.", "\n",
    "The genotype file contains", length(colnames(gt)), "individuals.", "\n",
    "The location file contains", nrow(locations), "individuals.", "\n")
cat("The filtering of phenotypic and genotypic file was done!", "\n", "\n")
#### Filter genotypes
filter_for_available_samples_in_gt <- function(i){
  if(any(locations$Sample_ID[i] == colnames(gt))){
    index <- which(locations$Sample_ID[i] == colnames(gt))
    return(index)}
}
index_gen <- mclapply(1:nrow(locations), filter_for_available_samples_in_gt, mc.cores = cores)
index_gen <- unlist(index_gen)
new_gen <- gt[, index_gen]
rownames(new_gen) <- SNP_IDs
gt <- new_gen
rm(new_gen)
cat("The phenotype file contains", length(colnames(phenotypes)), "individuals.", "\n",
    "The genotype file contains", length(colnames(gt)), "individuals.", "\n",
    "The location file contains", nrow(locations), "individuals.", "\n")
#### Filter phenotypes
colnames_gt <- colnames(gt)
filter_for_available_samples_in_phen <- function(i){
  if(any(colnames_gt[i] == colnames(phenotypes))){
    index <- which(colnames_gt[i] == colnames(phenotypes))
    return(index)}
}
index_phen <- mclapply(1:length(colnames_gt), filter_for_available_samples_in_phen, mc.cores = cores)
index_phen <- unlist(index_phen)
new_phen <- phenotypes[, index_phen]
phenotypes <- new_phen
#### Filter location data
filter_for_available_samples_in_location_data <- function(i){
  if(any(colnames(phenotypes)[i] == locations$Sample_ID)){
    index <- which(colnames(phenotypes)[i] == locations$Sample_ID)
    return(index)}
}
index_loc_data <- mclapply(1:length(phenotypes), filter_for_available_samples_in_location_data, mc.cores = cores)
index_loc_data <- unlist(index_loc_data)
new_loc <- locations[index_loc_data, ]
locations <- new_loc
ind_names_locations <- locations$Sample_ID
cat("The phenotype file contains", length(colnames(phenotypes)), "individuals.", "\n",
    "The genotype file contains", length(colnames(gt)), "individuals.", "\n",
    "The location file contains", nrow(locations), "individuals.", "\n")
index_phen <- mclapply(1:length(colnames_gt), filter_for_available_samples_in_phen, mc.cores = cores)
index_phen <- unlist(index_phen)
new_phen <- phenotypes[, index_phen]
phenotypes <- new_phen
cat("The phenotype file contains", length(colnames(phenotypes)), "individuals.", "\n",
    "The genotype file contains", length(colnames(gt)), "individuals.", "\n",
    "The location file contains", nrow(locations), "individuals.", "\n")
index_gen <- mclapply(1:nrow(locations), filter_for_available_samples_in_gt, mc.cores = cores)
index_gen <- unlist(index_gen)
new_gen <- gt[, index_gen]
rownames(new_gen) <- SNP_IDs
gt <- new_gen
rm(new_gen)
cat("The phenotype file contains", length(colnames(phenotypes)), "individuals.", "\n",
    "The genotype file contains", length(colnames(gt)), "individuals.", "\n",
    "The location file contains", nrow(locations), "individuals.", "\n")
cat("Names are:", colnames(phenotypes)[1:5], "\n",
    colnames(gt)[1:5],"\n",
    locations$Sample_ID[1:5],"\n")
index_loc_data <- mclapply(1:length(phenotypes), filter_for_available_samples_in_location_data, mc.cores = cores)
index_loc_data <- unlist(index_loc_data)
new_loc <- locations[index_loc_data, ]
locations <- new_loc
ind_names_locations <- locations$Sample_ID
cat("The phenotype file contains", length(colnames(phenotypes)), "individuals.", "\n",
    "The genotype file contains", length(colnames(gt)), "individuals.", "\n",
    "The location file contains", nrow(locations), "individuals.", "\n")
gt <- as.matrix(gt)
# LFMM -------------------------------------------------------------------------
cat("LFMM analysis starts.", "\n")
Y <- gt
Y <- t(as.matrix(Y))
latent_factor_K <- latent_factor_K
cat("Latent factor K was choosen to be:", "\t", latent_factor_K, "\n")
# Ridge estimates and GWAS test ------------------------------------------------
# Phenotype needs to be scaled
# Allele frequencies need to be imputed & 0 1 coded, phenotype and
# environmental variables needs to be scaled
### Phenotypes -----------------------------------------------------------------
phenotypes <- phenotypes[1:14,]
new_locations <- locations[,6:58]
new_locations <- t(new_locations)
rownames(new_locations) <- colnames(locations)[6:58]
colnames(new_locations) <- ind_names_locations
#phenotypes[1,which(is.na(phenotypes[1,]))] <- mean(as.numeric(phenotypes[1,]), na.rm = TRUE)
# Test function to see what is going wrong
### Environmental PCs ----------------------------------------------------------
new_env_PCs <- env_PCs[,6:58]
new_env_PCs <- t(new_env_PCs)
rownames(new_env_PCs) <- colnames(env_PCs)[6:58]
colnames(new_env_PCs) <- env_PCs$Sample_ID
cat("Calculate missingness at markers!", "\n")
calc_missingness <- function(dt){
  calc_missingness_per_marker <- function(i){
    nr_missing_obs <- length(which(Y[ ,i] == 9))
    rate_of_missing_obs <- nr_missing_obs/nrow(Y)
  }
  missingness_list <- mclapply(1:ncol(Y), calc_missingness_per_marker)
  missingness_dt <- unlist(missingness_list)
  missingness_dt <- matrix(data = missingness_dt, ncol = 1, nrow = ncol(Y))
  missingness_dt <- as.data.frame(missingness_dt)
  rownames(missingness_dt) <- colnames(Y)
  colnames(missingness_dt) <- "Rate_of_missingness"
  return(missingness_dt)
}
dt_missingness <- calc_missingness(dt = Y)
write.table(dt_missingness, paste0(result_dir, "Missingness_per_marker_before_filtering_", run_description, ".txt"))
cat("Remove markers with too high missing markers!", "\n")
missingness_value <- 0.05
remove_markers_to_high_missingness <- function(Y){
  replace_per_marker <- function(i){
    if(all(Y[ ,i] == 9)){
      index <- NA
      return(index)
    }
    if(length(which(Y[ ,i] == 9)) >= round((nrow(Y)*missingness_value), 0)){
      index <- NA
      return(index)
    }
    if(length(which(Y[ ,i] == 9)) <= round((nrow(Y)*missingness_value), 0)){
      index <- i
      return(index)
    }
  }
  index_mat <- mclapply(1:ncol(Y), replace_per_marker)
  index_mat <- unlist(index_mat)
  index <- which(!is.na(index_mat))
  Y <- Y[, index]
  return(Y)
}
new_Y <- remove_markers_to_high_missingness(Y = Y)
cat("Number of missing observations in the data set is:", length(which(new_Y==9)),"\n")
cat("Before we had", ncol(Y),"\n")
cat("Now we have", ncol(new_Y),"\n")
Y <- new_Y
rm(new_Y)
cat("Replace missing observations!", "\n")
# Start LFMM analysis ----------------------------------------------------------
Y = Y
#X_data = new_locations
X_data = new_env_PCs
latent_factor_K = latent_factor_K
cat("Following variable is processed:", "\t", rownames(X_data)[i], "\n")
var_name <- rownames(X_data)[i]
X <- X_data[i,]
X <- as.numeric(X)
X <- scale(X)
Y_only_var_sites <- Y[,which(colSums(Y)!=0)]
write.lfmm(Y_only_var_sites, paste0(lfmm_input_dir, var_name, "_genotypes.lfmm"))
write.lfmm(X, paste0(lfmm_input_dir, var_name, "_gradients.env"))
setwd(lfmm_input_dir)
#project.lfmm=lfmm(input.file = paste0(lfmm_input_dir, var_name, "_genotypes.lfmm"),
#                env = paste0(var_name, "_gradients.env"),
#                project = "new",
#                K = latent_factor_K,
#                iterations=100, burnin=50,
#                repetitions=10)
cat("Following variable is processed:", "\t", rownames(X_data)[i], "\n")
cat("Processing:", "\t", var_name, "\n")
cat("Project snmf starts.", "\n")
project.snmf <- snmf(input.file = paste0(lfmm_input_dir, var_name, "_genotypes.lfmm"),
                     K = latent_factor_K,
                     entropy = TRUE)
summary_snmf <- summary(project.snmf)
best = which.min(cross.entropy(project.snmf, K = latent_factor_K))
cat("Imputation starts.", "\n")
#impute(project.snmf, paste0(lfmm_input_dir, var_name, "_genotypes.lfmm"), method = 'mode', K = latent_factor_K, run = best)
impute(project.snmf, paste0(lfmm_input_dir, var_name, "_genotypes.lfmm"), method = 'mode', K = latent_factor_K, run = best)
cat("Replacing of missingness continues!", "\n")
Y <- read.lfmm(paste0(lfmm_input_dir, var_name, "_genotypes.lfmm_imputed.lfmm"))
#Y <- read.lfmm(paste0(lfmm_input_dir, "small_lfmm.lfmm"))
cat("Number of missing observations in the data set is:", length(which(Y==9)),"\n")
Y_mat <- matrix(data = Y,
                ncol = ncol(Y_only_var_sites),
                nrow = nrow(Y_only_var_sites),
                byrow = F)
#Y_mat <- matrix(data = Y,
#                ncol = 8,
#                nrow = nrow(Y_only_var_sites),
#                byrow = F)
cat("Matrix was built!", "\n")
print(Y_mat[1:4, 1:8])
rownames(Y_mat) <- rownames(Y_only_var_sites)
colnames(Y_mat) <- colnames(Y_only_var_sites)#[1:8]
replace_remaining_missing_obs <- function(Y){
  replace_per_marker <- function(i){
    gt_marker <- Y[ ,i]
    if(any(gt_marker == 9)){
      #cat("Missing observations found.", "\n")
      gt_marker[which(gt_marker == 9)] <- mean(as.numeric(gt_marker), na.rm = TRUE)
      new_gt <- as.numeric(gt_marker)
    }
    if(all(gt_marker != 9)){
      #cat("No missing observations.", "\n")
      new_gt <- as.numeric(gt_marker)
    }
    return(new_gt)
  }
  n_gt <- mclapply(1:ncol(Y), replace_per_marker, mc.cores = cores)
  n_gt <- unlist(n_gt)
  gt_new <- matrix(data = n_gt, nrow = nrow(Y), ncol = ncol(Y), byrow = FALSE)
  colnames(gt_new) <- colnames(Y)
  rownames(gt_new) <- rownames(Y)
  return(gt_new)
}
new_Y <- replace_remaining_missing_obs(Y = Y_mat)
cat("Remove non-variable sites!", "\n")
new_Y <- new_Y[,which(colSums(new_Y)!=0)]
cat("Run LFMM2!", "\n")
mod_lfmm2 <- lfmm2(input = new_Y,
                   env = X,
                   K = latent_factor_K,
                   effect.sizes = TRUE)
cat("LFMM2 was run!", "\n")
cat("P-values are getting calculated!", "\n")
pv <- lfmm2.test(object = mod_lfmm2,
                 input = new_Y,
                 env = X)
pvalues_pv <- pv$pvalues
effect_estimates <- mod_lfmm2@B
effect_estimates <- unlist(effect_estimates)
# Create QQ plot
png(filename = paste0(plot_dir, "QQplot_", Sys.Date(), "_",run_description, "_", rownames(X_data)[i], ".png"),
    width = 480, height = 480, units = "px")
qqplot(rexp(length(pvalues_pv), rate = log(10)),
       -log10(pvalues_pv), xlab = "Expected quantile",
       pch = 19, cex = .4)
qqline(rexp(length(pvalues_pv)))
dev.off()
# Create QQ plot based on adjusted value
png(filename = paste0(plot_dir, "QQplot_", Sys.Date(), "_adj_p_value", run_description, "_", rownames(X_data)[i], ".png"),
    width = 480, height = 480, units = "px")
qqplot(rexp(length(pvalues_pv), rate = log(10)),
       -log10(pvalues_pv), xlab = "Expected quantile",
       pch = 19, cex = .4)
qqline(rexp(length(pvalues_pv)))
dev.off()
# finished plotting
z_scores <- pv$zscores
dt_effect_estimates <- cbind(effect_estimates, pvalues_pv, z_scores)
new_dt_effect_estimates <- cbind(colnames(new_Y), rep(var_name, ncol(new_Y)), dt_effect_estimates)
ncol(new_dt_effect_estimates)
new_dt_effect_estimates <- as.data.frame(new_dt_effect_estimates)
print(head(new_dt_effect_estimates))
colnames(new_dt_effect_estimates) <- c("SNP_ID", "Variable", "effect_estimates", "p-values", "z-scores")
write.table(new_dt_effect_estimates, paste0(result_dir, Sys.Date(), "LFMM_based_on_", run_description, "_", var_name, ".txt"))
