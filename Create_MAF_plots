# MAF plots for PLINK GWAS results ---------------------------------------------
rm(list=ls())
cores <- 1
library(ggplot2)
library(foreach)
library(stringr)
library(parallel)
library(forcats)
## Load the data ---------------------------------------------------------------
setwd("/Beech_project/Romanian_Beech_data/GWAS/PLINK/MAF_plot/")
individual_names <- read.table("MAF_SNP_chr10_4898376_SD.nosex")
individual <- individual_names$V1
input_file_dir <- "/Beech_project/Romanian_Beech_data/Phenotypic data/New_phenotypic_data/"
setwd(input_file_dir)
locations <- read.table("2024_06_19_All_environmental_variables.txt", sep = "\t", header = TRUE)
pheno <- read.table("2024-08-30_RoBe_all_traits.txt", sep = "\t", header = TRUE)
plot_dir <- "/Beech_project/Romanian_Beech_data/GWAS/PLINK/MAF_plot/Plots/"
# Load plotting parameters -----------------------------------------------------
windowsFonts(my = windowsFont('Calibri'))
font_size <- 16
panel_background_colour <- "white"
plot_background_colour <- "white"
plot_font_colour <- "black"
colour_Lupului <- "maroon1"
colour_Ruia <- "skyblue1"
colour_Sinpetru <- "olivedrab4"
colour_Solomon <- "orange1"
colour_Tampa <- "indianred1"
colour_NA <- "grey45"
colour_populations <-c("Ruia" = colour_Ruia,
                       "Lupului" = colour_Lupului,
                       "Solomon" = colour_Solomon,
                       "Tampa" = colour_Tampa,
                       "Lempes" = colour_Sinpetru,
                       "NA" = colour_NA)
theme_paper <- theme(panel.background = element_rect(fill = panel_background_colour, colour = plot_background_colour),
                     panel.grid.minor = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.border =  element_rect(fill = NA, colour = panel_background_colour),
                     plot.background = element_rect(fill = plot_background_colour, colour = plot_background_colour),
                     axis.line = element_line(colour = plot_font_colour),
                     axis.line.x = element_line(color = plot_font_colour),
                     axis.ticks.y = element_blank(),
                     axis.ticks.x = element_line(color = plot_font_colour),
                     legend.position = "right",
                     legend.background =  element_rect(color = plot_font_colour, fill = panel_background_colour),
                     legend.box.background = element_rect(color = plot_font_colour, fill = panel_background_colour),
                     legend.title = element_text(colour = plot_font_colour, size = font_size, family = "my", face = "bold"),
                     legend.text = element_text(colour = plot_font_colour, size = font_size, family = "my"),
                     plot.title = element_text(colour = plot_font_colour, size = font_size, family = "my"),
                     axis.title.x = element_text(colour = plot_font_colour, size = font_size, family = "my", face = "bold"),
                     axis.text.x = element_text(size = font_size, family = "my", colour = plot_font_colour),
                     #axis.text.x = element_blank(),
                     axis.text.y = element_text(size = font_size, family = "my", colour = plot_font_colour),
                     axis.title.y = element_text(colour = plot_font_colour, size = font_size, family = "my", face = "bold"),
                     plot.tag = element_blank())
#panel_background_colour <- "white"
#plot_background_colour <- "white"
colour_populations <-c("Ruia" = colour_Ruia,
                       "Lupului" = colour_Lupului,
                       "Solomon" = colour_Solomon,
                       "Tampa" = colour_Tampa,
                       "Lempes" = colour_Sinpetru,
                       "NA" = colour_NA)
## Prepare the data ------------------------------------------------------------
index_SD <- which(colnames(pheno) == "Stomata_density")
phenotypes <- t(pheno[, index_SD])
colnames(phenotypes) <- pheno$Sample_ID
phenotypes <- as.data.frame(phenotypes)
plot_dir <- "/phD thesis/Presentation/Plots/"
## Read and plot MAF data ------------------------------------------------------
file_dir <- "/Beech_project/Romanian_Beech_data/GWAS/PLINK/MAF_plot/Results/"
prep_MAF_data <- function(file_dir){
  setwd(file_dir)
  file_names <- list.files()
  # Prepare the phenotypic data
  pheno <- t(phenotypes)
  pheno <- cbind(colnames(phenotypes), pheno)
  pheno <- as.data.frame(pheno)
  colnames(pheno) <- c("ID", "SD")
  pheno$Population <- 0
  pheno$Population <- gsub("[^a-zA-Z]", "", pheno$ID)
  # Calculate and trait mean
  mean_Si <- mean(as.numeric(pheno[which(pheno$Population == "Sinpetru"),2]), na.rm = TRUE)
  mean_Ta <- mean(as.numeric(pheno[which(pheno$Population == "Tampa"),2]), na.rm = TRUE)
  mean_So <- mean(as.numeric(pheno[which(pheno$Population == "Solomon"),2]), na.rm = TRUE)
  mean_Lu <- mean(as.numeric(pheno[which(pheno$Population == "Lupului"),2]), na.rm = TRUE)
  mean_Ru <- mean(as.numeric(pheno[which(pheno$Population == "Ruia"),2]), na.rm = TRUE)
  read_and_prep_file_per_file <- function(i){
    file_name <- file_names[i]
    marker_name <- str_sub(file_name, 9)
    new_marker_name <- gsub("._SD.raw", "", marker_name)
    cat("For marker", new_marker_name, "the MAF plot is getting prepared!", "\n")
    file <- read.table(file_name, header = TRUE)
    file <- as.data.frame(file)
    # Prepare the data
    filter_for_avail_samples <- function(i){
      if(any(colnames(phenotypes)[i] == file$IID)){
        index <- which(colnames(phenotypes)[i] == file$IID)
        return(index)}
    }
   index_data <- mclapply(1:length(colnames(phenotypes)), filter_for_avail_samples, mc.cores = cores)
   index_data <- unlist(index_data)
   new_file <- file[index_data, ]
   new_file$Populations <- 0
   new_file$Populations <- gsub("[^a-zA-Z]", "", new_file$IID)
   # Calculate the allele frequency
   Data_Si <- new_file[which(new_file$Populations == "Sinpetru"), 7]
   Data_Ta <- new_file[which(new_file$Populations == "Tampa"), 7]
   Data_So <- new_file[which(new_file$Populations == "Solomon"), 7]
   Data_Lu <- new_file[which(new_file$Populations == "Lupului"), 7]
   Data_Ru <- new_file[which(new_file$Populations == "Ruia"), 7]
   Af_Si <- sum(Data_Si, na.rm = TRUE)/(2*length(which(!is.na(Data_Si))))
   Af_Ta <- sum(Data_Ta, na.rm = TRUE)/(2*length(which(!is.na(Data_Ta))))
   Af_So <- sum(Data_So, na.rm = TRUE)/(2*length(which(!is.na(Data_So))))
   Af_Lu <- sum(Data_Lu, na.rm = TRUE)/(2*length(which(!is.na(Data_Lu))))
   Af_Ru <- sum(Data_Ru, na.rm = TRUE)/(2*length(which(!is.na(Data_Ru))))
   nr_of_total_ind_Si <- length(which(!is.na(Data_Si)))
   nr_of_total_ind_Ta <- length(which(!is.na(Data_Ta)))
   nr_of_total_ind_So <- length(which(!is.na(Data_So)))
   nr_of_total_ind_Lu <- length(which(!is.na(Data_Lu)))
   nr_of_total_ind_Ru <- length(which(!is.na(Data_Ru)))
   nr_HET_ind_Si <- length(which(Data_Si =="1"))
   nr_HET_ind_Ta <- length(which(Data_Ta =="1"))
   nr_HET_ind_So <- length(which(Data_So =="1"))
   nr_HET_ind_Lu <- length(which(Data_Lu =="1"))
   nr_HET_ind_Ru <- length(which(Data_Ru =="1"))
   nr_HOM_ind_Si <- length(which(Data_Si =="2"))
   nr_HOM_ind_Ta <- length(which(Data_Ta =="2"))
   nr_HOM_ind_So <- length(which(Data_So =="2"))
   nr_HOM_ind_Lu <- length(which(Data_Lu =="2"))
   nr_HOM_ind_Ru <- length(which(Data_Ru =="2"))
   dt_MAF_info <- rbind(cbind("Sinpetru", Af_Si, nr_HET_ind_Si, nr_HOM_ind_Si, nr_of_total_ind_Si),
                        cbind("Tampa", Af_Ta, nr_HET_ind_Ta, nr_HOM_ind_Ta, nr_of_total_ind_Ta),
                        cbind("Solomon", Af_So, nr_HET_ind_So, nr_HOM_ind_So, nr_of_total_ind_So),
                        cbind("Lupului", Af_Lu, nr_HET_ind_Lu, nr_HOM_ind_Lu, nr_of_total_ind_Lu),
                        cbind("Ruia", Af_Ru, nr_HET_ind_Ru, nr_HOM_ind_Ru, nr_of_total_ind_Ru))
   dt_MAF_info <- as.data.frame(dt_MAF_info)
   colnames(dt_MAF_info) <- c("Stand", "MAF","Nr_HET", "Nr_HOM", "Nr_Tot")
   dt_MAF_info$MAF <- round(as.numeric(dt_MAF_info$MAF), 4)
   cat("Allele freq: in Sinpetru", Af_Si, "\n",
       "Allele freq: in Tampa", Af_Ta, "\n",
       "Allele freq: in Solomon", Af_So, "\n",
       "Allele freq: in Lupului", Af_Lu, "\n",
       "Allele freq: in Ruia", Af_Ru, "\n")
   dt_MAF <- cbind(rep(new_marker_name, 5),
                   c("Lempes", "Tampa", "Solomon", "Lupului", "Ruia"),
                   rbind(Af_Si, Af_Ta, Af_So, Af_Lu, Af_Ru),
                   rbind(mean_Si, mean_Ta, mean_So, mean_Lu, mean_Ru))
   dt_MAF <- as.data.frame(dt_MAF)
   colnames(dt_MAF) <- c("Marker", "Population", "Allele_freq", "Mean_SD")
   dt_MAF$Mean_SD <- as.numeric(dt_MAF$Mean_SD)
   dt_MAF$Allele_freq <- as.numeric(dt_MAF$Allele_freq)
   cor_coeff <- cor(dt_MAF$Allele_freq, dt_MAF$Mean_SD)
   cor_test <- cor.test(dt_MAF$Allele_freq, dt_MAF$Mean_SD)#
   p_val <- cor_test$p.value
   dt_info_cor <- cbind(new_marker_name, dt_MAF_info, cor_coeff, p_val)
   # Plotting
   MAF_plot <- ggplot(data = dt_MAF, aes(y = Allele_freq, x = Mean_SD))+
     theme_paper +
     geom_smooth(method = lm, se = TRUE, colour = plot_font_colour)+
     geom_point(aes(y = Allele_freq, x = Mean_SD,
                    colour = fct_inorder(Population)), size = 2)+
     labs(x = "Stomatal density", y = "MAF", tag = "A")+
     scale_colour_manual(values = colour_populations, name = "Stand")
   # Save plot
   ggsave(paste0(plot_dir, Sys.Date(), "MAF_plot_", new_marker_name,".png"),
          MAF_plot,
          device = png,
          height = 3,
          width = 5,
          dpi = 900)
   return(dt_info_cor)
  }
 dt_combined <- foreach(i = 1:length(file_names), .combine = rbind) %do% read_and_prep_file_per_file(i)
 return(dt_combined)
}
dt_MAF_combined <- prep_MAF_data(file_dir = "/Beech_project/Romanian_Beech_data/GWAS/PLINK/MAF_plot/Results/")
write.table(dt_MAF_combined , paste0("/Beech_project/Romanian_Beech_data/GWAS/PLINK/MAF_plot/MAF_all_sig_markers_for_SD.txt"))
dt_MAF_combined <- as.data.frame(dt_MAF_combined)
range(as.numeric(dt_MAF_combined$p_val))
round(range(as.numeric(dt_MAF_combined$cor_coeff)),4)
# Plot GT against phenotyp -----------------------------------------------------
theme_paper <- theme(panel.background = element_rect(fill = panel_background_colour, colour = plot_background_colour),
                     panel.grid.minor = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.border =  element_rect(fill = NA, colour = panel_background_colour),
                     plot.background = element_rect(fill = plot_background_colour, colour = plot_background_colour),
                     axis.line = element_line(colour = plot_font_colour),
                     axis.line.x = element_line(color = plot_font_colour),
                     axis.ticks.y = element_blank(),
                     axis.ticks.x = element_blank(),
                     legend.position = "none",
                     legend.box.background = element_rect(color = plot_font_colour),
                     legend.title = element_text(colour = plot_font_colour, size = font_size, family = "my", face = "bold"),
                     legend.text = element_text(colour = plot_font_colour, size = font_size, family = "my"),
                     plot.title = element_text(colour = plot_font_colour, size = font_size, family = "my"),
                     axis.title.x = element_text(colour = plot_font_colour, size = font_size, family = "my", face = "bold"),
                     axis.text.x = element_blank(),
                     #axis.text.x = element_blank(),
                     axis.text.y = element_text(size = font_size, family = "my", colour = plot_font_colour),
                     axis.title.y = element_text(colour = plot_font_colour, size = font_size, family = "my", face = "bold"),
                     plot.tag = element_blank())
colour_populations <-c("Ruia" = colour_Ruia,
                       "Lupului" = colour_Lupului,
                       "Solomon" = colour_Solomon,
                       "Tampa" = colour_Tampa,
                       "Sinpetru" = colour_Sinpetru,
                       "NA" = colour_NA)
prep_GT_data_and_plot <- function(file_dir){
  plot_dir <- "/Beech_project/Romanian_Beech_data/GWAS/PLINK/GT_plots/"
  setwd("/Beech_project/Romanian_Beech_data/GWAS/PLINK/MAF_plot/Results/")
  file_names <- list.files()
  # Prepare the phenotypic data
  pheno <- t(phenotypes)
  pheno <- cbind(colnames(phenotypes), pheno)
  pheno <- as.data.frame(pheno)
  colnames(pheno) <- c("ID", "SD")
  pheno$Population <- 0
  pheno$Population <- gsub("[^a-zA-Z]", "", pheno$ID)
  get_GT_plot_per_marker <- function(i){
    file_name <- file_names[i]
    marker_name <- str_sub(file_name, 9)
    new_marker_name <- gsub("._SD.raw", "", marker_name)
    cat("For marker", new_marker_name, "the MAF plot is getting prepared!", "\n")
    file <- read.table(file_name, header = TRUE)
    file <- as.data.frame(file)
    # Prepare the data
    filter_for_avail_samples <- function(i){
      if(any(colnames(phenotypes)[i] == file$IID)){
        index <- which(colnames(phenotypes)[i] == file$IID)
        return(index)}
    }
    index_data <- mclapply(1:length(colnames(phenotypes)), filter_for_avail_samples, mc.cores = cores)
    index_data <- unlist(index_data)
    new_file <- file[index_data, ]
    new_file$Populations <- 0
    new_file$Populations <- gsub("[^a-zA-Z]", "", new_file$IID)
    filter_for_avail_samples <- function(i){
      if(any(file$IID[i] == colnames(phenotypes))){
        index <- which(file$IID[i] == colnames(phenotypes))
        return(index)}
    }
    index_data <- mclapply(1:length(file$IID), filter_for_avail_samples, mc.cores = cores)
    index_data <- unlist(index_data)
    pheno <- pheno[index_data, ]
    # Calculate the allele frequency
    Data_Si <- new_file[which(new_file$Populations == "Sinpetru"), ]
    Data_Ta <- new_file[which(new_file$Populations == "Tampa"), ]
    Data_So <- new_file[which(new_file$Populations == "Solomon"), ]
    Data_Lu <- new_file[which(new_file$Populations == "Lupului"), ]
    Data_Ru <- new_file[which(new_file$Populations == "Ruia"), ]
    dt <- rbind(Data_Si, Data_Ta, Data_So, Data_Lu, Data_Ru)
    pheno <- pheno[order(pheno$ID), ]
    dt <- dt[order(dt$FID), ]
    dt <- cbind(dt, pheno)
    dt <- dt[, c(2,7,8,10)]
    colnames(dt) <- c("Sample_ID", "Genotype", "Population", "Stomatal_density")
    dt$Genotype <- as.numeric(dt$Genotype)
    dt$Stomatal_density <- as.numeric(dt$Stomatal_density)
    dt <- dt[which(!is.na(dt$Genotype)),]
    lm_regression <- lm(formula = Stomatal_density ~ Genotype, data = dt)
    sum_regr_coeff <- summary(lm_regression)
    reg_coeff <- sum_regr_coeff$r.squared
    cat("Marker:", new_marker_name, "\n")
    print(sum_regr_coeff)
    # Plotting
    GT_plot <- ggplot(data = dt, aes(y = Stomatal_density, x = Genotype))+
      theme_paper +
      geom_smooth(method = lm, se = TRUE, colour = "black")+
      geom_jitter(aes(y = Stomatal_density, x = Genotype,
                    colour = fct_inorder(Population)), size = 1, width = 0.1)+
      labs(x = "Genotypes", y = "Stomatal density")+
      scale_colour_manual(values = colour_populations, name = "Stand")#+
    scale_x_discrete(labels = c("0/0", "HET", "1/1"))
    ## Save plot
    ggsave(paste0(plot_dir, Sys.Date(), "_GT_plot_", new_marker_name,".png"),
         GT_plot,
         device = png,
         height = 4,
         width = 4,
         dpi = 900)
    dt <- cbind(new_marker_name, reg_coeff)
    dt <- as.data.frame(dt)
    return(dt)
  }
  dt_combined <- foreach(i = 1:length(file_names), .combine = rbind) %do% get_GT_plot_per_marker(i)
  return(dt_combined)
}
dt_GT_combined <- prep_GT_data_and_plot(file_dir = "/Beech_project/Romanian_Beech_data/GWAS/PLINK/MAF_plot/Results/")
dt_GT_combined <- as.data.frame(dt_GT_combined)
dt_GT_combined$reg_coeff <- as.numeric(dt_GT_combined$reg_coeff)
dt_GT_combined$reg_coeff <- round(dt_GT_combined$reg_coeff, 4)
dt_GT_combined[which(dt_GT_combined$reg_coeff > 0.23),]
write.table(dt_GT_combined, paste0("/Beech_project/Romanian_Beech_data/GWAS/PLINK/MAF_plot/Results/Regression_and_correlation_coefficients.txt"))
