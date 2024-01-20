### The BEGINNING ~~~~~
##
# ~ Plots NLSparrow--ngsAdmix | First written by Jose Samaniego with later modifications by George Pacheco.


# Cleans the environment ~ 
rm(list=ls())


# Sets working directory ~
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Loads required packages ~
pacman::p_load(tidyverse, scales, optparse, plyr, RColorBrewer, extrafont, gtable, grid, mdthemes, ggtext, glue)


# Creates colour palette ~
nb.cols <- 15
MyColours <- colorRampPalette(brewer.pal(8, "Paired"))(nb.cols)


# Loads the data ~
samples <- read.table("NLSparrow--ngsAdmix.popfile", stringsAsFactors = FALSE, sep = "\t")


# Reads the annotation file ~
ids <- read.table("NLSparrow.labels", stringsAsFactors = FALSE, sep = "\t", header = FALSE)


# Adds column ids names ~
colnames(ids) <- c("Sample_ID")


# Expands ids by adding Population ~
ids$Population <- ifelse(grepl("FR0", ids$Sample_ID), "Sales",
                  ifelse(grepl("KAZ", ids$Sample_ID), "Chokpak",
                  ifelse(grepl("Lesina", ids$Sample_ID), "Lesina",
                  ifelse(grepl("Crotone", ids$Sample_ID), "Crotone",
                  ifelse(grepl("Guglionesi", ids$Sample_ID), "Guglionesi",
                  ifelse(grepl("PI22NLD0001M", ids$Sample_ID), "Target",
                  ifelse(grepl("PD22NLD0146F", ids$Sample_ID), "Garderen",
                  ifelse(grepl("PD22NLD0147F", ids$Sample_ID), "Garderen",
                  ifelse(grepl("PDOM2022NLD0", ids$Sample_ID), "Utrecht", "Error")))))))))


# Reorders Population ~
ids$Population <- factor(ids$Population, ordered = T,
                         levels = c("Utrecht",
                                    "Garderen",
                                    "Sales",
                                    "Crotone",
                                    "Guglionesi",
                                    "Lesina",
                                    "Chokpak",
                                    "Target"))


# Expands PCA_Annot by adding Species ~
ids$Species <- ifelse(ids$Population %in% c("Utrecht", "Sales", "Garderen"), "House",
               ifelse(ids$Population %in% c("Chokpak", "Lesina"), "Spanish",
               ifelse(ids$Population %in% c("Guglionesi", "Crotone"), "Italian",
               ifelse(ids$Population %in% c("Target"), "Uncertain", "Error"))))


# Reorders Population ~
ids$Species <- factor(ids$Species, ordered = T,
                      levels = c("House",
                                 "Italian",
                                 "Spanish",
                                 "Uncertain"))


# Ask Sama ~
fulldf <- data.frame()


# Ask Sama 2 ~
x <- list(c(4, 5, 7, 6, 3, 2, 1),
          c(1, 3, 5, 6, 2, 4),
          c(5, 3, 2, 1, 4),
          c(2, 4, 3, 1),
          c(1, 2 ,3),
          c(1, 2))


# Defines samples' IDs ~
sampleid = "Sample_ID"


# Ask Sama 3 ~
for (j in 1:length(samples[,1])){
  data <- read.table(samples[j,1])[, x[[j]]]
  for (i in 1:dim(data)[2]) { 
    temp <- data.frame(Ancestry = data[, i])
    temp$K <- as.factor(rep(i, times = length(temp$Ancestry)))
    temp[sampleid] <- as.factor(ids[sampleid][,1])
    temp$K_Value <- as.factor(rep(paste("K = ", dim(data)[2], sep = ""), times = length(temp$Ancestry)))
    temp <- merge(ids, temp)
    fulldf <- rbind(fulldf, temp)}}


# Defines the target to be plotted ~
target = "Population"


# Creates the plots ~
ngsAdmix <-
 ggplot(fulldf, aes(x = Sample_ID, y = Ancestry, fill = K)) +
 geom_bar(stat = "identity", width = .85, alpha = .70) +
 facet_grid(K_Value ~ get(target), space = "free_x", scales = "free_x") +
 scale_fill_manual(values = c("#ee0000", "#1E90FF", "#FFD700", "#88419d", "#9ecae1", "#c6dbef", "#fec44f")) +
 scale_x_discrete(expand = c(0, 0)) + 
 scale_y_continuous(expand = c(0, 0), breaks = NULL) +
 theme(panel.background = element_rect(fill = "#ffffff"),
       panel.grid.minor.x = element_blank(),
       panel.grid.major = element_blank(),
       panel.spacing = unit(.3, "lines"),
       plot.title = element_blank(),
       legend.position = "none",
       axis.title = element_blank(),
       axis.text.x.bottom = element_text(colour = "#000000", face = "bold", angle = 90, vjust = .5, hjust = .5),
       #axis.text.x.bottom = element_blank(),
       axis.text.y = element_blank(),
       axis.ticks = element_blank(),
       strip.background = element_rect(colour = "#000000", fill = "#FAFAFA", linewidth = .05),
       strip.text.x = element_text(colour = "#000000", face = "bold", size = 16, angle = 90, margin = margin(.5, .1, .5, .1, "cm")),
       strip.text.y = element_text(colour = "#000000", face = "bold", size = 16, angle = 90, margin = margin(0, .1, 0, .1, "cm")))
 
 
 # Saves the final plot ~
 ggsave(ngsAdmix, file = "NLSparrow.MinInd90.Autosomes_TMP_SNPs_XXX.pdf",
        device = cairo_pdf, width = 16, height = 8, dpi = 600)


# Adds grob ~
ngsAdmix_G <- ggplotGrob(ngsAdmix)
ngsAdmix_G <- gtable_add_rows(ngsAdmix_G, unit(1.25, "cm"), pos = 5)


# Adds top strips ~
ngsAdmix_G <- gtable_add_grob(ngsAdmix_G, list(rectGrob(gp = gpar(col = "#000000", fill = "#1E90FF", alpha = .7, size = .75, lwd = .25)),
               textGrob("HOUSE", gp = gpar(cex = 1.4, fontface = 'bold', col = "black"))),
               t = 6, l = 4, b = 6, r = 9, name = c("a", "b"))
ngsAdmix_G <- gtable_add_grob(ngsAdmix_G, list(rectGrob(gp = gpar(col = "#000000", fill = "#FFD700", alpha = .7, size = .5, lwd = .25)),
               textGrob("ITALIAN", gp = gpar(cex = 1.4, fontface = 'bold', col = "black"))),
               t = 6, l = 11, b = 6, r = 13, name = c("a", "b"))
ngsAdmix_G <- gtable_add_grob(ngsAdmix_G, list(rectGrob(gp = gpar(col = "#000000", fill = "#ee0000", alpha = .7, size = .75, lwd = .25)),
               textGrob("SPANISH", gp = gpar(cex = 1.4, fontface = 'bold', col = "black"))),
               t = 6, l = 15, b = 6, r = 17, name = c("a", "b"))


# Controls separation ~
ngsAdmix_G <- gtable_add_rows(ngsAdmix_G, unit(2/10, "line"), 6)


# Creates the final plot ~
grid.newpage()
grid.draw(ngsAdmix_G)


# Saves the final plot ~
ggsave(ngsAdmix_G, file = "AllSamples_haplotypecaller.raw.vcf.Filtered.MAF20.OnlyAutosomes.pdf",
       device = cairo_pdf, width = 28, height = 12, scale = .9, dpi = 600)
ggsave(ngsAdmix_G, file = "AllSamples_haplotypecaller.raw.vcf.Filtered.MAF20.OnlyAutosomes.jpeg",
       width = 28, height = 12, scale = .9, dpi = 600)


#
##
### The END ~~~~~



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20),
#c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19),
#c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18),
#c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),
#c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16),
#c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15),
#c(1,2,3,4,5,6,7,8,9,10,11,12,13,14),
#c(1,2,3,4,5,6,7,8,9,10,11,12,13),
#c(1,2,3,4,5,6,7,8,9,10,11,12),
#c(1,2,3,4,5,6,7,8,9,10,11),
#c(1,2,3,4,5,6,7,8,9,10),
#c(1,2,3,4,5,6,7,8,9),
#c(1,2,3,4,5,6,7,8),
#c(1,2,3,4,5,6,7),
#c(1,2,3,4,5,6),
#c(1,2,3,4,5),
#c(1,2,3,4),
#c(1,2,3),
#c(1,2))