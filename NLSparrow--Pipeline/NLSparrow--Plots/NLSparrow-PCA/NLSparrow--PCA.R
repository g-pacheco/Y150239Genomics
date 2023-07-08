### The BEGINNING ~~~~~
##
# ~ Plots NLSparrow--PCA | First written by Hom?re J. Alves Monteiro with later modifications by George Pacheco.


# Cleans the environment ~ 
rm(list=ls())


# Sets working directory ~
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Loads required packages ~
pacman::p_load(optparse, tidyverse, plyr, RColorBrewer, extrafont, ggforce, ggstar, ggrepel, RcppCNPy, reshape2, lemon, plotly,
               gridExtra, grid, cowplot, patchwork, ggpubr, rphylopic)


# Imports extra fonts ~
font_import()


# Loads data ~
data <- as.matrix(read.table("NLSparrow_MinInd90_SNPs.cov"), stringsAsFactors = FALSE)
annot <- read.table("NLSparrow.labels", sep = "\t", header = FALSE, stringsAsFactors = FALSE)


# Runs PCA ~
PCA <- eigen(data)


# Merges the first 3 PCs with annot ~
PCA_Annot <- as.data.frame(cbind(annot, PCA$vectors[, c(1:3)]))
colnames(PCA_Annot) <- c("Sample_ID", "PCA_1", "PCA_2", "PCA_3")

# Expands PCA_Annot by adding Population ~
PCA_Annot$Population <- #ifelse(grepl("8L", PCA_Annot$Sample_ID), "Helgeland",
                        #ifelse(grepl("8M", PCA_Annot$Sample_ID), "Helgeland",
                        ifelse(grepl("FR0", PCA_Annot$Sample_ID), "Sales",
                        ifelse(grepl("KAZ", PCA_Annot$Sample_ID), "Chokpak",
                        ifelse(grepl("Lesina", PCA_Annot$Sample_ID), "Lesina",
                        ifelse(grepl("Guglionesi", PCA_Annot$Sample_ID), "Guglionesi",
                        ifelse(grepl("PI22NLD0001M", PCA_Annot$Sample_ID), NA,
                        ifelse(grepl("PD22NLD0146F", PCA_Annot$Sample_ID), "Garderen",
                        ifelse(grepl("PD22NLD0147F", PCA_Annot$Sample_ID), "Garderen",
                        ifelse(grepl("PDOM2022NLD0", PCA_Annot$Sample_ID), "Utrecht", "Error"))))))))


# Reorders Population ~
PCA_Annot$Population <- factor(PCA_Annot$Population, ordered = T,
                        levels = c(#"Helgeland",
                                   "Utrecht",
                                   "Sales",
                                   "Garderen",
                                   "Guglionesi",
                                   "Lesina",
                                   "Chokpak",
                                   NA))


# Expands PCA_Annot by adding Species ~
PCA_Annot$Species <- ifelse(PCA_Annot$Population %in% c("Utrecht", "Sales", "Garderen"), "House",
                     ifelse(PCA_Annot$Population %in% c("Chokpak", "Lesina"), "Spanish",
                     ifelse(PCA_Annot$Population %in% c("Guglionesi"), "Italian",
                     ifelse(PCA_Annot$Population %in% NA, NA, "Error"))))


# Reorders Population ~
PCA_Annot$Species <- factor(PCA_Annot$Species, ordered = T,
                               levels = c("House",
                                          "Italian",
                                          "Spanish",
                                          NA))


# Gets Eigenvalues of each Eigenvectors ~
PCA_Eigenval_Sum <- sum(PCA$values)
(PCA$values[1]/PCA_Eigenval_Sum)*100
(PCA$values[2]/PCA_Eigenval_Sum)*100
(PCA$values[3]/PCA_Eigenval_Sum)*100


# Defines the shapes to be used for each Group ~
Shapes <- as.vector(c(1, 2, 3, 13, 11, 23))


# Creates legend plot ~
MyLegend_Plot <-
  ggplot(data = PCA_Annot, aes_string(x = "PCA_1", y = "PCA_2")) +
  geom_star(aes(starshape = Population, fill = Species), size = 2.8, starstroke = .15, alpha = .85) +
  scale_fill_manual(values = c("#1E90FF", "#FFD700", "#ee0000"), na.translate = FALSE) +
  scale_starshape_manual(values = Shapes, na.translate = F) +
  scale_x_continuous("PC 1 (4.3%)",
                     #breaks = c(0.99, 1, 1.01),
                     #labels = c("0.99", "1", "1.01"),
                     #limits = c(-.37, .145),
                     expand = c(.005, .005)) +
  scale_y_continuous("PC 2 (1.2%)",
                     #breaks = c(-.15, 0, .15, .3, .45), 
                     #labels = c("-.15", "0", ".15", ".3", ".45"), 
                     #limits = c(-0.0525, 0.0525),
                     expand = c(.03, .03)) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        legend.position = "top",
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.margin = margin(t = 0, b = 0, r = 15, l = 15),
        legend.box = "vertical",
        legend.box.margin = margin(t = 20, b = 30, r = 0, l = 0)) +
  guides(starshape = guide_legend(title = "Population", title.theme = element_text(size = 16, face = "bold"),
                                  label.theme = element_text(size = 15),
                                  override.aes = list(starshape = Shapes, size = 5, starstroke = .15), nrow = 1, order = 2),
         fill = guide_legend(title = "Species", title.theme = element_text(size = 16, face = "bold"),
                             label.theme = element_text(size = 15),
                             override.aes = list(starshape = 21, size = 5, starstroke = .15), nrow = 1, order = 1),
         colour = "none")


# Saves plot ~
#ggsave(MyLegend_Plot, file = "XXX.pdf",
#       device = cairo_pdf, limitsize = FALSE, scale = .65, width = 16, height = 16, dpi = 600)


# Defines the shapes to be used for each Group ~
Shapes_2 <- as.vector(c(1, 2, 3, 13, 11, 23, 14))


# Combines all populations from the Faroe Islands ~
PCA_Annot$Species <- as.character(PCA_Annot$Species)
PCA_Annot$Population <- as.character(PCA_Annot$Population)
PCA_Annot <- PCA_Annot %>%
  mutate_at(c("Population", "Species"), ~replace_na(., "Target"))


#geom_label_repel(data = Coords_Global, size = 3.75, seed = 10, min.segment.length = 0, force = 25, segment.curvature = 1,
#                 nudge_x = 0, nudge_y = 0, max.overlaps = Inf, fontface = "bold", colour = "black",
#                 aes(x = Longitude, y = Latitude, label = LocationOnly,
#                 fill = Class_Article), alpha = 0.9, show.legend = FALSE) +


# Reorders Population ~
PCA_Annot$Population <- factor(PCA_Annot$Population, ordered = T,
                               levels = c(#"Helgeland",
                                          "Utrecht",
                                          "Sales",
                                          "Garderen",
                                          "Guglionesi",
                                          "Lesina",
                                          "Chokpak",
                                          "Target"))


# Reorders Population ~
PCA_Annot$Species <- factor(PCA_Annot$Species, ordered = T,
                            levels = c("House",
                                       "Italian",
                                       "Spanish",
                                       "Target"))


#PCA_AnnotLabel <-  PCA_Annot %>%
#  filter((Species == "Target"))

PCA_AnnotAbridged <-  PCA_Annot %>%
  filter((Sample_ID != "PDOM2022NLD0046M" & Sample_ID != "PDOM2022NLD0083F"))


# Expands PCA_Annot by adding Labels ~
PCA_Annot$Labels <- ifelse(PCA_Annot$Species %in% c("Target"), "Target\nIndividual", "")


PCA_12_Blog <-
  ggplot(data = PCA_AnnotAbridged, aes_string(x = "PCA_1", y = "PCA_2")) +
  geom_star(aes(starshape = Population, fill = Species), alpha = .85, size = 2.8, starstroke = .15) +
  scale_fill_manual(values = c("#1E90FF", "#FFD700", "#ee0000", "#d9d9d9")) +
  scale_starshape_manual(values = Shapes_2) +
  geom_label_repel(data = PCA_Annot, aes(label = Labels),
                   family = ".SF Compact Rounded", size = 3.8, fontface = "bold", max.overlaps = 100, nudge_x = .055, nudge_y = .05,
                   point.padding = .6, segment.size = .3, colour = "black", fill = "#d9d9d9", alpha = .85,
                   arrow = arrow(angle = 30, length = unit(.10, "inches"),
                   ends = "last", type = "open")) +
  geom_mark_ellipse(aes(filter = Species == "House", label = "House\nSparrow"), con.colour = "#1E90FF", colour = "#1E90FF",
                    label.fill = "#d9d9d9", expand = unit(4, "mm"), con.border = "one", label.fontsize = 10.65,
                    con.type = "straight", label.family = ".SF Compact Rounded", con.cap = 0, label.hjust = .5, show.legend = FALSE) +
  geom_mark_ellipse(aes(filter = Species == "Spanish", label = "Spanish\nSparrow"), con.colour = "#ee0000", colour = "#ee0000",
                    label.fill = "#d9d9d9", expand = unit(4, "mm"), con.border = "one", label.fontsize = 10.65,
                    con.type = "elbow", label.family = ".SF Compact Rounded", con.cap = 0, label.hjust = .5, show.legend = FALSE) +
  geom_mark_ellipse(aes(filter = Species == "Italian", label = "Italian\nSparrow"), con.colour = "#FFD700", colour = "#FFD700",
                    label.fill = "#d9d9d9", expand = unit(4, "mm"), con.border = "one", label.fontsize = 10.65,
                    con.type = "elbow", label.family = ".SF Compact Rounded", con.cap = 0, label.hjust = .5, show.legend = FALSE) +
  scale_x_continuous("PC 1 (10.80%)",
                     #breaks = c(0.99, 1, 1.01),
                     #labels = c("0.99", "1", "1.01"),
                     #limits = c(-.10, .17),
                     expand = c(.02, .02)) +
  scale_y_continuous("PC 2 (2.75%)",
                     #breaks = c(-.08, -.04, 0.00), 
                     #labels = c("-0.08", "-0.04", "0.00"),
                     limits = c(-.11, .03),
                     expand = c(.02, .02)) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        legend.position = "none",
        axis.title.x = element_text(size = 15, face = "bold", margin = margin(t = 18, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 15, face = "bold", margin = margin(t = 0, r = 18, b = 0, l = 0)),
        axis.text = element_text(color = "#000000", size = 12, face = "bold"),
        axis.ticks = element_line(color = "#000000", linewidth = .3),
        strip.text = element_text(colour = "#000000", size = 13, face = "bold"),
        strip.background = element_rect(colour = "#000000", fill = "#d6d6d6", linewidth = .3),
        axis.line = element_line(colour = "#000000", linewidth = .3))


#ggsave(PCA_12_Blog, file = "YYY2.pdf",
#       device = cairo_pdf, limitsize = FALSE, scale = .65, width = 16, height = 16, dpi = 600)


# Isolates legend ~
MyLegendBlog <- get_legend(MyLegend_Plot)


# Gets final plot ~
PCA_Plot_Blog <- ggarrange(PCA_12_Blog, nrow = 1, legend.grob = MyLegendBlog)


# Saves plot ~
ggsave(PCA_Plot_Blog, file = "NLSparrow_Blog_Abridged-Ellipses--PCA.pdf",
       device = cairo_pdf, limitsize = FALSE, scale = 1, width = 9, height = 9, dpi = 600)


PCA_12 <-
  ggplot(data = PCA_Annot, aes(x = "PCA_1", y = "PCA_2")) +
  geom_star(aes(starshape = Species, fill = Population), size = 2.8, starstroke = .15) +
  scale_fill_manual(values = c("#fee8c8", "#fdbb84", "#a6bddb", "#2b8cbe", "#31a354", "#dd3497", "#fa9fb5")) +
  scale_starshape_manual(values = Shapes) +
  scale_x_continuous("PC 1 (13.82%)",
                     #breaks = c(0.99, 1, 1.01),
                     #labels = c("0.99", "1", "1.01"),
                     #limits = c(-.3625, .14),
                     expand = c(.005, .005)) +
  scale_y_continuous("PC 2 (3.06%)",
                     #breaks = c(-.2, 0, .2, .4), 
                     #labels = c("-.2", "0", ".2", ".4"),
                     #limits = c(-0.0525, 0.0525),
                     expand = c(.03, .03)) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        legend.position = "none",
        axis.title.x = element_text(size = 15, face = "bold", margin = margin(t = 15, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 15, face = "bold", margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.text = element_text(color = "#000000", size = 12, face = "bold"),
        axis.ticks = element_line(color = "#000000", linewidth = .3),
        strip.text = element_text(colour = "#000000", size = 13, face = "bold"),
        strip.background = element_rect(colour = "#000000", fill = "#d6d6d6", linewidth = .3),
        axis.line = element_line(colour = "#000000", linewidth = .3))


PCA_13 <-
  ggplot(data = PCA_Annot, aes_string(x = "PCA_1", y = "PCA_3")) +
  geom_star(aes(starshape = Species, fill = Population), size = 2.8, starstroke = .15) +
  scale_fill_manual(values = c("#fee8c8", "#fdbb84", "#a6bddb", "#2b8cbe", "#31a354", "#dd3497", "#fa9fb5")) +
  scale_starshape_manual(values = Shapes) +
  scale_x_continuous("PC 1 (13.82)",
                     #breaks = c(0.99, 1, 1.01),
                     #labels = c("0.99", "1", "1.01"),
                     #limits = c(-.195, .106),
                     expand = c(.005, .005)) +
  scale_y_continuous("PC 3 (2.45%)",
                     #position = "right",
                     #breaks = c(-0.05, -0.025, 0, 0.025, 0.05), 
                     #labels = c("-0.05", "-0.025", "0", "0.025", "0.05"), 
                     #limits = c(-0.0525, 0.0525),
                     expand = c(.03, .03)) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        legend.position = "none",
        axis.title.x = element_text(size = 15, face = "bold", margin = margin(t = 15, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 15, face = "bold", margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.text = element_text(color = "#000000", size = 12, face = "bold"),
        axis.ticks = element_line(color = "#000000", linewidth = .3),
        strip.text = element_text(colour = "#000000", size = 13, face = "bold"),
        strip.background = element_rect(colour = "#000000", fill = "#d6d6d6", linewidth = .3),
        axis.line = element_line(colour = "#000000", linewidth = .3))


PCA_23 <-
  ggplot(data = PCA_Annot, aes_string(x = "PCA_2", y = "PCA_3")) +
  geom_star(aes(starshape = Species, fill = Population), size = 2.8, starstroke = .15) +
  scale_fill_manual(values = c("#fee8c8", "#fdbb84", "#a6bddb", "#2b8cbe", "#31a354", "#dd3497", "#fa9fb5")) +
  scale_starshape_manual(values = Shapes) +
  scale_x_continuous("PC 2 (3.06%)",
                     #breaks = c(0.99, 1, 1.01),
                     #labels = c("0.99", "1", "1.01"),
                     #limits = c(-.25, .1275),
                     expand = c(.005, .005)) +
  scale_y_continuous("PC 3 (2.45%)",
                     #position = "right",
                     #breaks = c(-0.05, -0.025, 0, 0.025, 0.05), 
                     #labels = c("-0.05", "-0.025", "0", "0.025", "0.05"), 
                     #limits = c(-0.0525, 0.0525),
                     expand = c(.03, .03)) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        legend.position = "none",
        axis.title.x = element_text(size = 15, face = "bold", margin = margin(t = 15, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 15, face = "bold", margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.text = element_text(color = "#000000", size = 12, face = "bold"),
        axis.ticks = element_line(color = "#000000", linewidth = .3),
        strip.text = element_text(colour = "#000000", size = 13, face = "bold"),
        strip.background = element_rect(colour = "#000000", fill = "#d6d6d6", linewidth = .3),
        axis.line = element_line(colour = "#000000", linewidth = .3))


# Isolates legend ~
MyLegend <- get_legend(MyLegend_Plot)


# Gets final plot ~
PCA_Plot <- ggarrange(PCA_12, PCA_13, PCA_23, nrow = 3, legend.grob = MyLegend)


# Saves plot ~
ggsave(PCA_Plot, file = "NLSparrow--PCA.pdf",
       device = cairo_pdf, limitsize = FALSE, scale = .75, width = 16, height = 16, dpi = 600)


# Creates plot ~
PCA_Plot <-
ggplot(PCA_Annot, aes_string(x = "PCA_1", y = "PCA_2")) +
  geom_star(aes(starshape = Species, fill = Population), size = 2.75, colour = "#000000", alpha = .9) +
  scale_fill_manual(values = c("#fef0d9", "#fdd49e", "#fdbb84", "#fc8d59", "#e34a33", "#b30000", "#990000")) +
  scale_x_continuous("PC 1 (13.82%)",
                     #breaks = c(0.99, 1, 1.01),
                     #labels = c("0.99", "1", "1.01"),
                     #limits = c(-0.275, 0.15),
                     expand = c(0, 0)) +
  scale_y_continuous("PC 2 (3.06%)",
                     #breaks = c(-0.05, -0.025, 0, 0.025, 0.05), 
                     #labels = c("-0.05", "-0.025", "0", "0.025", "0.05"), 
                     #limits = c(-0.0525, 0.0525),
                     expand = c(.015, .015)) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.position = "right",
        legend.title = element_text(color = "#000000", size = 13),
        legend.text = element_text(size = 11),
        axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 18, face = "bold", margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text = element_text(color = "#000000", size = 13),
        axis.ticks = element_line(color = "#000000", size = .3),
        axis.line = element_line(colour = "#000000", linewidth = .3)) +
  guides(fill = guide_legend(title = "Population", title.theme = element_text(size = 15, face = "bold"),
                             label.theme = element_text(size = 14)),
         fill = guide_legend(title = "Population", title.theme = element_text(size = 15, face = "bold"),
                             label.theme = element_text(size = 14)))


# Saves plot ~
ggsave(PCA_Plot, file = "NLSparrow--PCA.pdf",
       device = cairo_pdf, scale = 1, width = 12, height = 8, dpi = 600)


#
##
### The END ~~~~~