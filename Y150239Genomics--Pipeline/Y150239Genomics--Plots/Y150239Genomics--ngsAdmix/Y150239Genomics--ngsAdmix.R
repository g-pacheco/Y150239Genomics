### The BEGINNING ~~~~~
##
# ~ Plots Y150239Genomics--ngsAdmix | Written by Jose Samaniego & George Pacheco.


# Cleans the environment ~ 
rm(list=ls())


# Sets working directory ~
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Loads required packages ~
pacman::p_load(tidyverse, scales, optparse, plyr, RColorBrewer, extrafont, gtable, grid, mdthemes, ggtext, glue, ggh4x)


# Creates colour palette ~
nb.cols <- 15
MyColours <- colorRampPalette(brewer.pal(8, "Paired"))(nb.cols)


# Loads the data ~
samples.auto <- read.table("Y150239Genomics--ngsAdmix.Autosomes.popfile", stringsAsFactors = FALSE, sep = "\t")
samples.allo <- read.table("Y150239Genomics--ngsAdmix.Allosome.popfile", stringsAsFactors = FALSE, sep = "\t")


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


# Expands PCA_Annot by adding Species ~
ids$Species <- ifelse(ids$Population %in% c("Utrecht", "Sales", "Garderen", "Meerkerk"), "House",
               ifelse(ids$Population %in% c("Chokpak", "Lesina"), "Spanish",
               ifelse(ids$Population %in% c("Guglionesi", "Crotone"), "Italian",
               ifelse(ids$Population %in% c("Y150239"), "Y150239", "Error"))))


# Reorders Population ~
ids$Species <- factor(ids$Species, ordered = T,
                      levels = c("House",
                                 "Italian",
                                 "Spanish",
                                 "Y150239"))


# Recognises chromosomes ~
ids.auto <- ids
ids.allo <- ids
ids.auto$chrtype <- "Autosomes"
ids.allo$chrtype <- "Allosome"


# Creates data frame ~
fulldf.auto <- data.frame()
fulldf.allo <- data.frame()


x.auto <- list(c(3, 2, 1, 4, 6, 5, 7),
               c(4, 5, 1, 2, 6, 3),
               c(4, 2, 5, 3, 1),
               c(4, 3, 1, 2),
               c(1, 3, 2),
               c(1, 2))

x.allo <- list(c(5, 6, 3, 1, 4, 2, 7),
               c(6, 2, 5, 4, 3, 1),
               c(3, 4, 2, 5, 1),
               c(3, 4, 2, 1),
               c(3, 1, 2),
               c(1, 2))


# Defines samples' IDs ~
sampleid = "Sample_ID"


# Loops over all Ks while adding labels and reordering clusters ~
for (j in 1:length(samples.auto[, 1])){
  data <- read.table(samples.auto[j, 1])[, x.auto[[j]]]
  for (i in 1:dim(data)[2]) { 
    temp <- data.frame(Ancestry = data[, i])
    temp$K <- as.factor(rep(i, times = length(temp$Ancestry)))
    temp[sampleid] <- as.factor(ids.auto[sampleid][, 1])
    temp$K_Value <- as.factor(rep(paste("K = ", dim(data)[2], sep = ""), times = length(temp$Ancestry)))
    temp <- merge(ids.auto, temp)
    fulldf.auto <- rbind(fulldf.auto, temp)}}


for (j in 1:length(samples.allo[, 1])){
  data <- read.table(samples.allo[j, 1])[, x.allo[[j]]]
  for (i in 1:dim(data)[2]) { 
    temp <- data.frame(Ancestry = data[, i])
    temp$K <- as.factor(rep(i, times = length(temp$Ancestry)))
    temp[sampleid] <- as.factor(ids.allo[sampleid][, 1])
    temp$K_Value <- as.factor(rep(paste("K = ", dim(data)[2], sep = ""), times = length(temp$Ancestry)))
    temp <- merge(ids.allo, temp)
    fulldf.allo <- rbind(fulldf.allo, temp)}}


fulldf <- rbind(fulldf.auto, fulldf.allo)


# Reorders chrtye ~
fulldf$chrtype <- factor(fulldf$chrtype, ordered = T,
                        levels = c("Autosomes",
                                    "Allosome"))


# Defines the target to be plotted ~
target = "Population"


# Creates the plot ~
ngsAdmix <-
ggplot(fulldf, aes(x = Sample_ID, y = Ancestry, fill = K)) +
  geom_bar(stat = "identity", width = .85, alpha = .70) +
  facet_nested(chrtype + K_Value ~ Population, scales = "free_x", space = "free_x",
               strip = strip_nested(background_x = list(element_rect(colour = "#000000", fill = "#FAFAFA", linewidth = .05)),
                                    background_y = list(element_rect(colour = "#000000", fill = "#bdbdbd", linewidth = .05),
                                                        element_rect(colour = "#000000", fill = "#FAFAFA", linewidth = .05)), by_layer_y = TRUE)) +
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
        #axis.text.x.bottom = element_text(colour = "#000000", face = "bold", angle = 90, vjust = .5, hjust = .5),
        axis.text.x.bottom = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        strip.text.x = element_text(colour = "#000000", face = "bold", family = "Optima", size = 14, angle = 90, margin = margin(.5, .1, .5, .1, "cm")),
        strip.text.y = element_text(colour = "#000000", face = "bold", family = "Optima", size = 14, angle = 90, margin = margin(0, .1, 0, .1, "cm")))


# Saves the final plot ~
ggsave(ngsAdmix, file = "XXX.pdf",
       device = cairo_pdf, width = 20, height = 12, dpi = 600)


# Adds grob ~
ngsAdmix_G <- ggplotGrob(ngsAdmix)
ngsAdmix_G <- gtable_add_rows(ngsAdmix_G, unit(1.25, "cm"), pos = 5)


# Adds top strips ~
ngsAdmix_G <- gtable_add_grob(ngsAdmix_G, list(rectGrob(gp = gpar(col = "#000000", fill = "#1E90FF", alpha = .7, size = .75, lwd = .25)),
               textGrob("HOUSE", gp = gpar(cex = 1.5, fontface = 'bold', fontfamily = "Optima", col = "black"))),
               t = 6, l = 4, b = 6, r = 9, name = c("a", "b"))
ngsAdmix_G <- gtable_add_grob(ngsAdmix_G, list(rectGrob(gp = gpar(col = "#000000", fill = "#FFD700", alpha = .7, size = .5, lwd = .25)),
               textGrob("ITALIAN", gp = gpar(cex = 1.5, fontface = 'bold', fontfamily = "Optima", col = "black"))),
               t = 6, l = 11, b = 6, r = 13, name = c("a", "b"))
ngsAdmix_G <- gtable_add_grob(ngsAdmix_G, list(rectGrob(gp = gpar(col = "#000000", fill = "#ee0000", alpha = .7, size = .75, lwd = .25)),
               textGrob("SPANISH", gp = gpar(cex = 1.5, fontface = 'bold', fontfamily = "Optima", col = "black"))),
               t = 6, l = 15, b = 6, r = 17, name = c("a", "b"))


# Controls separation ~
ngsAdmix_G <- gtable_add_rows(ngsAdmix_G, unit(2 / 10, "line"), 6)


# Creates the final plot ~
grid.newpage()
grid.draw(ngsAdmix_G)


# Saves the final plot ~
ggsave(ngsAdmix_G, file = "Leiden_Admix.pdf",
       device = cairo_pdf, width = 30, height = 12, scale = 1, dpi = 600)
ggsave(ngsAdmix_G, file = "Leiden_Admix.jpeg",
       width = 30, height = 12, scale = 1, dpi = 600)


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