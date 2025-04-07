# Plot GWAS PLINK results ------------------------------------------------------
rm(list = ls())
library(data.table)
library(ggplot2)
#install.packages("forcats")
library(forcats)
library(stringr)
library(dplyr)
library(parallel)
### Test data on local machine -------------------------------------------------
setwd("/Beech_project/Romanian_Beech_data/GWAS/PLINK/GWAS_results/")
# SD
GWAS_results = fread("SD_MAFMindfilt.assoc.linear", head=T)
GWAS_perm = fread("SD_MAFMindfilt_perm.assoc.linear.perm", head=T)
var_name <- "SD"
# D13 2020
var_name <- "D13_2020"
GWAS_results = fread("D13_2020_MAFMindfilt_weather_data_2020_2021_2022_noSept.assoc.linear", head=T)
GWAS_perm = fread("D13_2020_MAFMindfilt_weather_data_2020_2021_2022_noSept_perm.assoc.linear.perm", head=T)
# D13 2021
var_name <- "D13_2021"
GWAS_results = fread("D13_2021_MAFMindfilt_weather_data_2020_2021_2022_noSept.assoc.linear", head=T)
GWAS_perm = fread("D13_2021_MAFMindfilt_weather_data_2020_2021_2022_noSept_perm.assoc.linear.perm", head=T)
# D13 2022
var_name <- "D13_2022"
GWAS_results = fread("D13_2022_MAFMindfilt_weather_data_2020_2021_2022_noSept.assoc.linear", head=T)
GWAS_perm = fread("D13_2022_MAFMindfilt_weather_data_2020_2021_2022_noSept_perm.assoc.linear.perm", head=T)
### Prepare the data -----------------------------------------------------------
GWAS_results <- as.data.frame(GWAS_results)
GWAS_results$SNP <- paste0(GWAS_results$CHR, "_",GWAS_results$BP)
length(which(GWAS_results$P < 0.01))
length(which(GWAS_results$P < 0.001))
length(which(GWAS_results$P < 0.00001))
length(which(GWAS_results$P < 0.00001))
length(which(GWAS_results$P < 0.000001))

sig_results <- GWAS_results[(which(GWAS_results$P < 0.000001)),]
SNP_markers <- unique(sig_results$SNP)
tab <- matrix(data = unlist(str_split(SNP_markers, "_")), byrow = T, ncol = 2)
table(tab[,1])

length(unique(GWAS_results[which(GWAS_results$P < 0.00001),2]))
GWAS_results[(which(GWAS_results$P < 0.00001 & GWAS_results$CHR == 10)),]
# Permuted results
length(which(GWAS_perm$EMP1 < 0.001))
length(which(GWAS_perm$EMP1 < 0.0001))
length(which(GWAS_perm$EMP1 < 0.00001))
length(which(GWAS_perm$EMP1 < 0.000005))
length(which(GWAS_perm$EMP1 < 0.000001))
length(which(GWAS_results$P < 0.000001))
# Significance threshold
sig_threshold <- 0.000001
# Save significant markers -----------------------------------------------------
sig_marker <- GWAS_results[which(GWAS_results$P < sig_threshold),]
table(sig_marker$CHR)
length(unique(sig_marker$SNP))
write.table(sig_marker, paste0("C:/Users/mtost/Documents/Beech_project/Romanian_Beech_data/GWAS/PLINK/Significant markers/", Sys.Date(),"NEW_2025_Significant_markers_", var_name, "_", sig_threshold, ".txt"))
### Prepare data set -----------------------------------------------------------
GWAS_results$CHR <- as.numeric(GWAS_results$CHR)
GWAS_results$SNP <- paste0(GWAS_results$CHR, "_", GWAS_results$BP)
GWAS_results$SNP <- as.factor(GWAS_results$SNP)
head(GWAS_results)
GWAS_results <- GWAS_results[order(GWAS_results$CHR), ]
GWAS_results$CHR <- as.factor(GWAS_results$CHR)
GWAS_results$log_p <- -log10(GWAS_results$P)
GWAS_results <- GWAS_results[which(GWAS_results$SNP != "10_10907267" & GWAS_results$TEST != "DOMDEV"),]
### Load plotting parameters ---------------------------------------------------
windowsFonts(my = windowsFont('Calibri'))
font_size <- 14
plot_background_colour <- "white"
plot_font_colour <- "black"
# For presentation
plot_background_colour <- "#00003e"
plot_font_colour <- "white"
font_size <- 16
plot_background_colour <- "#00003e"
plot_font_colour <- "white"
theme_Manhattan <- theme(panel.background = element_rect(fill = plot_background_colour, colour = plot_background_colour),
                         panel.grid.minor = element_blank(),
                         panel.grid.major = element_blank(),
                         panel.border =  element_rect(fill = NA, colour = plot_background_colour),
                         plot.background = element_rect(fill = plot_background_colour, colour = plot_background_colour),
                         axis.line = element_line(colour = plot_font_colour),
                         axis.line.x = element_line(color = plot_font_colour),
                         axis.ticks.x = element_blank(),
                         axis.ticks.y = element_line(color = plot_font_colour),
                         legend.position = "bottom",
                         legend.text = element_text(size = font_size, family = "my", colour = plot_font_colour),
                         legend.title =  element_text(colour =plot_font_colour, size = font_size, family = "my", face = "bold"),
                         plot.title = element_text(hjust = 0.5, colour = plot_font_colour, size = font_size),
                         axis.title.x = element_text(hjust = 0.5, colour =plot_font_colour,family = "my", size = font_size, face = "bold"),
                         axis.text.x = element_blank(),
                         axis.text.y = element_text(size = font_size, family = "my", colour = plot_font_colour),
                         axis.title.y = element_text(hjust = 0.5, colour = plot_font_colour, size = font_size, family = "my", face = "bold"),
                         plot.tag = element_text(hjust = 0.5, colour =plot_font_colour,family = "my", size = font_size))

my_pal_col_1 <- (c("1" =  "#00003e", "2" = "chocolate4",
                   "3" = "#6699CC", "4" = "darkorchid4",
                   "5" = "#F0E442", "6" = "darkcyan",
                   "7" = "#661100", "8" = "darkseagreen3",
                   "9" = "#CC79A7", "10" = "darkorange",
                   "11" = "royalblue2", "12" = "#882255"))
my_pal_col_1 <- (c("1" =  "aliceblue", "2" = "chocolate4",
                 "3" = "#6699CC", "4" = "darkorchid4",
                 "5" = "#F0E442", "6" = "darkcyan",
                 "7" = "#661100", "8" = "darkseagreen3",
                 "9" = "#CC79A7", "10" = "darkorange",
                 "11" = "royalblue2", "12" = "#882255"))
sig_threshold <- 0.000001
# Manhattan plots --------------------------------------------------------------
colnames(GWAS_results)
GWAS_results$CHR <- as.character(GWAS_results$CHR)
GWAS_results$SNP
GWAS_results[which(GWAS_results$CHR == "10" & GWAS_results$BP > 13661222),]
plot_Manhattan <- function(data,
                           my_pal_col){
  Manhattan_plot <- ggplot(data = data,
                           aes(x = fct_inorder(SNP),
                               y = log_p,
                               colour = fct_inorder(CHR)))+
    geom_point()+
    theme_Manhattan +
    labs(x = "Position (bp)", y = "-log10(p-value)")+
    geom_hline(yintercept = -log10(sig_threshold), colour = "red")+
    scale_color_manual(values = my_pal_col, name = "Chromosomes")+
    geom_vline(xintercept = "10_4896732", colour = "red", linewidth = 0.5)+
    geom_vline(xintercept = "10_13677756", colour = "red", linewidth = 0.5)+
    ylim(0, 10)
  return(Manhattan_plot)
}
Manhattan_plot_A <- plot_Manhattan(data = GWAS_results,
                                   my_pal_col = my_pal_col_1)
ggsave(paste0(plot_dir, Sys.Date(),"_MAFMissFilt_Manhattanplot_", var_name,".png"),
       Manhattan_plot_A,
       device = png,
       height = 6,
       width = 10,
       dpi = 600)
### QQ-Plot --------------------------------------------------------------------
# Prepare the data
GWAS_results$CHR <- as.numeric(GWAS_results$CHR)
GWAS_results$SNP <- paste0(GWAS_results$CHR, "_", GWAS_results$BP)
GWAS_results$SNP <- as.factor(GWAS_results$SNP)
head(GWAS_results)
GWAS_results <- GWAS_results[order(GWAS_results$CHR), ]
GWAS_results$CHR <- as.factor(GWAS_results$CHR)
GWAS_results$log_p <- -log10(GWAS_results$P)
### Prepare the data -----------------------------------------------------------
GWAS_results <- as.data.frame(GWAS_results)
GWAS_results$SNP <- paste0(GWAS_results$CHR, "_",GWAS_results$BP)
theme_QQPlot <- theme(panel.background = element_rect(fill = plot_background_colour, colour = plot_background_colour),
                      panel.grid.minor = element_blank(),
                      panel.grid.major = element_blank(),
                      panel.border =  element_rect(fill = NA, colour = plot_background_colour),
                      plot.background = element_rect(fill = plot_background_colour, colour = plot_background_colour),
                      axis.line = element_line(colour = plot_font_colour),
                      axis.line.x = element_line(color = plot_font_colour),
                      axis.ticks.x = element_line(color = plot_font_colour),
                      axis.ticks.y = element_line(color = plot_font_colour),
                      legend.position = "none",
                      legend.text = element_text(size = font_size, family = "my", colour = plot_font_colour),
                      legend.title =  element_text(colour =plot_font_colour, size = font_size, family = "my", face = "bold"),
                      plot.title = element_text(hjust = 0.5, colour = plot_font_colour, size = font_size),
                      axis.title.x = element_text(hjust = 0.5, colour =plot_font_colour,family = "my", size = font_size, face = "bold"),
                      axis.text.x = element_text(size = font_size, family = "my", colour = plot_font_colour),
                      axis.text.y = element_text(size = font_size, family = "my", colour = plot_font_colour),
                      axis.title.y = element_text(hjust = 0.5, colour = plot_font_colour, size = font_size, family = "my", face = "bold"),
                      plot.tag = element_text(hjust = 0.5, colour =plot_font_colour,family = "my", size = font_size))
plot_QQ_plot <- function(data,
                         Plot_nr){
  QQ_plot <- ggplot(data = data,
                    aes(sample = log_p))+
    geom_qq(colour = plot_font_colour)+
    theme_QQPlot+
    labs(tag = Plot_nr, x = "Theoretical quantiles", y = "Sample quantiles")+
    #geom_qq_line()#+
    #stat_qq_line(colour = "steelblue", fullrange = T)+
    geom_abline(slope=1,intercept = 0, colour = "steelblue")+
    ylim(0,5)+
    xlim(0,5)
  return(QQ_plot)
}
# SD
QQ_plot <- plot_QQ_plot(data = GWAS_results,
                        Plot_nr = "a")
# NC2020
QQ_plot <- plot_QQ_plot(data = GWAS_results,
                        Plot_nr = "b")
# NC2021
QQ_plot <- plot_QQ_plot(data = GWAS_results,
                        Plot_nr = "c")
# NC2022
QQ_plot <- plot_QQ_plot(data = GWAS_results,
                        Plot_nr = "d")
# CNR 2020
QQ_plot <- plot_QQ_plot(data = GWAS_results,
                        Plot_nr = "e")
# CNR 2021
QQ_plot <- plot_QQ_plot(data = GWAS_results,
                        Plot_nr = "f")
# CNR 2022
QQ_plot <- plot_QQ_plot(data = GWAS_results,
                        Plot_nr = "g")
# D13C 2020
QQ_plot <- plot_QQ_plot(data = GWAS_results,
                        Plot_nr = "h")
# D13C 2021
QQ_plot <- plot_QQ_plot(data = GWAS_results,
                        Plot_nr = "i")
# D13C 2022
QQ_plot <- plot_QQ_plot(data = GWAS_results,
                        Plot_nr = "j")
### Save plots -----------------------------------------------------------------
ggsave(paste0("C:/Users/mtost/Documents/Beech_project/Romanian_Beech_data/Write_up/Plots/PLINK_GWAS_Manhattan/NEW/", Sys.Date(),"_PLINK_QQplot_", var_name,".png"),
       QQ_plot,
       device = png,
       height = 6,
       width = 6,
       dpi = 600)
