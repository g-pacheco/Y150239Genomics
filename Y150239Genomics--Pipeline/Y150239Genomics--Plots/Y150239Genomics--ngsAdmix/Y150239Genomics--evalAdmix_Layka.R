### The BEGINNING ~~~~~
##
# ~ Plots NLSparrow--evalAdmix | First written by Jose Samaniego with later modifications by George Pacheco.


# Cleans the environment ~ 
rm(list=ls())


# Sets working directory ~
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Loads required packages ~
pacman::p_load(tidyverse, pheatmap, scales, optparse, plyr, RColorBrewer, extrafont, gtable, grid, mdthemes, ggtext, glue)
source("visFuns.R")


# Loads the data ~
Corres_2 <- as.matrix(read.table("AllSamples_bcftools.raw.vcf.Filtered.Allosome.NoKinship.NoTreeSparrow.MAFfiltered.Pruned.K2.corres")) 
Corres_3 <- as.data.frame(read.table("AllSamples_bcftools.raw.vcf.Filtered.Allosome.NoKinship.NoTreeSparrow.MAFfiltered.Pruned.K3.corres"))
Corres_4 <- as.data.frame(read.table("AllSamples_bcftools.raw.vcf.Filtered.Allosome.NoKinship.NoTreeSparrow.MAFfiltered.Pruned.K4.corres"))
Corres_5 <- as.data.frame(read.table("AllSamples_bcftools.raw.vcf.Filtered.Allosome.NoKinship.NoTreeSparrow.MAFfiltered.Pruned.K5.corres"))
Corres_6 <- as.data.frame(read.table("AllSamples_bcftools.raw.vcf.Filtered.Allosome.NoKinship.NoTreeSparrow.MAFfiltered.Pruned.K6.corres"))
Corres_7 <- as.data.frame(read.table("AllSamples_bcftools.raw.vcf.Filtered.Allosome.NoKinship.NoTreeSparrow.MAFfiltered.Pruned.K7.corres"))


# Reads the annotation file ~
ids <- read.table("Y150239Genomics--ngsAdmix.labels", stringsAsFactors = FALSE, sep = "\t", header = FALSE)


# Adds column ids names ~
colnames(ids) <- c("Sample_ID")


# Expands ids by adding Population ~
ids$Population <- ifelse(grepl("FR0", ids$Sample_ID), "Sales",
                  ifelse(grepl("KAZ", ids$Sample_ID), "Chokpak",
                  ifelse(grepl("Lesina", ids$Sample_ID), "Lesina",
                  ifelse(grepl("Crotone", ids$Sample_ID), "Crotone",
                  ifelse(grepl("Guglionesi", ids$Sample_ID), "Guglionesi",
                  ifelse(grepl("PI22NLD0001M", ids$Sample_ID), "Y150239",
                  ifelse(grepl("PD22NLD0146F", ids$Sample_ID), "Garderen",
                  ifelse(grepl("PD22NLD0147F", ids$Sample_ID), "Garderen",
                  ifelse(grepl("PDOM2022NLD0077M", ids$Sample_ID), "Meerkerk",
                  ifelse(grepl("PDOM2022NLD0", ids$Sample_ID), "Utrecht", "Error"))))))))))


# Reorders Population ~
ids$Population <- factor(ids$Population, ordered = T,
                         levels = c("Utrecht",
                                    "Garderen",
                                    "Meerkerk",
                                    "Sales",
                                    "Crotone",
                                    "Guglionesi",
                                    "Lesina",
                                    "Chokpak",
                                    "Y150239"))


# Plot correlation of residuals ~
plotCorRes(cor_mat = Corres_2, pop = as.vector(ids[, 2]), ord = NULL,
           title = "Evaluation of 1000G admixture proportions with K = 2",
           plot_legend = T, max_z = 0.1, min_z = -.1)





Kinship_Plot_Heatmap <-
  ggplot(Corres_2, aes(Ind1, Ind2, fill = rab)) + 
  geom_tile(colour = "#000000") +
  scale_fill_continuous(low = "#ffffff", high = "#f768a1") +
  #scale_fill_viridis(trans = "reverse", option = "plasma") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(limits = rev, expand = c(0, 0)) +
  facet_wrap(Population ~., scales = "free", ncol = 2) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(1, "lines"),
        legend.position = "right",
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.margin = margin(t = 0, b = 0, r = 15, l = 15),
        legend.box = "vertical",
        legend.box.margin = margin(t = 20, b = 30, r = 0, l = 0),
        axis.title = element_blank(),
        axis.text.x = element_text(color = "#000000", size = 8.25, face = "bold", angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(color = "#000000", size = 8.25, face = "bold"),
        axis.ticks = element_line(color = "#000000", linewidth = .3),
        strip.text = element_text(colour = "#000000", size = 14, face = "bold", family = "Optima"),
        strip.background = element_rect(colour = "#000000", fill = "#d6d6d6", linewidth = .3),
        axis.line = element_line(colour = "#000000", linewidth = .3)) +
  guides(fill = guide_legend(title = "Rab", title.theme = element_text(size = 16, face = "bold"),
                             label.theme = element_text(size = 15), reverse = TRUE))




plotCorRes(cor_mat = Corres_3, pop = as.vector(ids[, 2]), ord = NULL,
           title = "Evaluation of 1000G admixture proportions with K = 3", plot_legend = T, max_z = 0.1, min_z = -0.1)

plotCorRes(cor_mat = Corres_4, pop = as.vector(ids[, 2]), ord = NULL,
           title = "Evaluation of 1000G admixture proportions with K = 4", plot_legend = T, max_z = 0.1, min_z = -0.1)

plotCorRes(cor_mat = Corres_5, pop = as.vector(ids[, 2]), ord = NULL,
           title = "Evaluation of 1000G admixture proportions with K = 5", plot_legend = T, max_z = 0.1, min_z = -0.1)

plotCorRes(cor_mat = Corres_6, pop = as.vector(ids[, 2]), ord = NULL,
           title = "Evaluation of 1000G admixture proportions with K = 6", plot_legend = T, max_z = 0.1, min_z = -0.1)

plotCorRes(cor_mat = Corres_7, pop = as.vector(ids[, 2]), ord = NULL,
           title = "Evaluation of 1000G admixture proportions with K = 7", plot_legend = T, max_z = 0.1, min_z = -0.1)


# Creates & Saves the heatmap ~
pheatmap(Corres_2, cluster_rows = TRUE, cluster_cols = TRUE, border_color = "black", cellwidth = 10, cellheight = 10,
         filename = "K2_T.pdf")
pheatmap(Corres_3, cluster_rows = FALSE, cluster_cols = FALSE, border_color = "black", cellwidth = 10, cellheight = 10,
         filename = "K3.pdf")
pheatmap(Corres_4, cluster_rows = FALSE, cluster_cols = FALSE, border_color = "black", cellwidth = 10, cellheight = 10,
         filename = "K4.pdf")
pheatmap(Corres_5, cluster_rows = FALSE, cluster_cols = FALSE, border_color = "black", cellwidth = 10, cellheight = 10,
         filename = "K5.pdf")
pheatmap(Corres_6, cluster_rows = FALSE, cluster_cols = FALSE, border_color = "black", cellwidth = 10, cellheight = 10,
         filename = "K6.pdf")
pheatmap(Corres_7, cluster_rows = FALSE, cluster_cols = FALSE, border_color = "black", cellwidth = 10, cellheight = 10,
         filename = "K7.pdf")



#
##
### The END ~~~~~