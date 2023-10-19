### The BEGINNING ~~~~~
##
# ~ Plots FPG--MDS | First written by Filipe G. Vieira with later modifications by George Pacheco.


# Cleans the environment ~ 
rm(list=ls())


# Sets working directory ~
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Loads required packages ~
pacman::p_load(optparse, tidyverse, plyr, RColorBrewer, extrafont, ggforce, ggstar, ggrepel, RcppCNPy, reshape2, lemon, plotly,
               gridExtra, grid, cowplot, patchwork, ggpubr, rphylopic)


# Imports extra fonts ~
#loadfonts(device = "win", quiet = TRUE)


# Loads list of options ~
option_list <- list(make_option(c('-i','--in_file'), action = 'store', type = 'character', default = "stdin", help = 'Input file'),
                    make_option(c('-a','--annot'), action = 'store', type = 'character', default = NULL, help = 'File with indiv annotations'),
                    make_option(c('--id_column'), action = 'store', type = 'numeric', default = 1, help = 'Column to use as ID'),
                    make_option(c('-L','--in_maj_labels'), action = 'store', type = 'character', default = NULL, help = 'Column from annotation file to use as MAJOR label'),
                    make_option(c('-l','--in_min_labels'), action = 'store', type = 'character', default = NULL, help = 'Column from annotation file to use as MINOR label'),
                    make_option(c('--no_header'), action = 'store_true', type = 'logical', default = FALSE, help = 'Input file has no header'),
                    make_option(c('--var_excl'), action = 'store', type = 'character', default = NULL, help = 'Variables to exclude from analysis'))

  
# Defines parameters ~
opt <- parse_args(OptionParser(option_list = option_list))
opt$in_file = "AllSamples_haplotypecaller.raw.vcf.Filtered.vcf2beagle.mds"
opt$annot = "NLSparrow.labels"
opt$id_column = 1
opt$in_maj_labels = "Population"


# Reads data ~
data <- read.table(opt$in_file, row.names = 1, sep = "\t", header = !opt$no_header, stringsAsFactors = FALSE, check.names = FALSE)
n <- ncol(data)


# Reads annotation file ~
if(!is.null(opt$annot)){
annot <- read.table(opt$annot, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
colnames(annot) <- c("Sample_ID")


# Expands PCA_Annot by adding Population ~
annot$Population <- ifelse(grepl("FR0", annot$Sample_ID), "Sales",
                    ifelse(grepl("KAZ", annot$Sample_ID), "Chokpak",
                    ifelse(grepl("Lesina", annot$Sample_ID), "Lesina",
                    ifelse(grepl("Crotone", annot$Sample_ID), "Crotone",
                    ifelse(grepl("Guglionesi", annot$Sample_ID), "Guglionesi",
                    ifelse(grepl("PI22NLD0001M", annot$Sample_ID), NA,
                    ifelse(grepl("PD22NLD0146F", annot$Sample_ID), "Garderen",
                    ifelse(grepl("PD22NLD0147F", annot$Sample_ID), "Garderen",
                    ifelse(grepl("PDOM2022NLD0", annot$Sample_ID), "Utrecht", "Error")))))))))


# Reorders Population ~
annot$Population <- factor(annot$Population, ordered = T,
                          levels = c("Utrecht",
                                     "Sales",
                                     "Garderen",
                                     "Crotone",
                                     "Guglionesi",
                                     "Lesina",
                                     "Chokpak",
                                     NA))



data <- merge(data, annot, by.x = 0, by.y = opt$id_column)


# Gets rownames back into place ~
rownames(data) <- data[,1]; data <- data[,-1]
data[colnames(annot)[opt$id_column]] <- rownames(data)}


# Excludes variables ~
if( !is.null(opt$var_excl)){
  opt$var_excl <- unlist(strsplit(opt$var_excl, ","))
  data <- data[!(rownames(data) %in% opt$var_excl),]}


# Gets Major labels mean location ~
colors <- NULL
if(!is.null(opt$in_maj_labels)){
  
  
  # Merge Major labels ~
  in_maj_labels <- unlist(strsplit(opt$in_maj_labels, ",", fixed = TRUE))
  tmp_data <- data[,in_maj_labels[1]]
  data[in_maj_labels[1]] <- NULL
  if(length(in_maj_labels) > 1){
    for (cnt in 2:length(in_maj_labels)){
      tmp_data <- paste(tmp_data, data[,in_maj_labels[cnt]], sep="/")
      data[in_maj_labels[cnt]] <- NULL}
    opt$in_maj_labels <- "MERGE"}
  
  
  # Make sure Major label column is after data ~
  data <- data.frame(data, tmp_data)
  colnames(data)[ncol(data)] <- opt$in_maj_labels
  
  
  # Converts to factor, in case there is a Major label with just numbers~
  data[,opt$in_maj_labels] <- factor(data[,opt$in_maj_labels])
  
  
  # If label was in input file, decreases number of data columns ~
  if(is.null(opt$annot) || !opt$in_maj_labels %in% colnames(annot))
    n = n - 1
  
  
  # Gets mean value for Major label ~
  data_mean <- ddply(data, opt$in_maj_labels, function(x){colMeans(x[, 1:n], na.rm = TRUE)})
  colors <- as.character(opt$in_maj_labels)}


# Expands PCA_Annot by adding Species ~
data$Species <- ifelse(data$Population %in% c("Utrecht", "Sales", "Garderen"), "House",
                ifelse(data$Population %in% c("Chokpak", "Lesina"), "Spanish",
                ifelse(data$Population %in% c("Crotone", "Guglionesi"), "Italian",
                ifelse(data$Population %in% NA, NA, "Error"))))


# Reorders Population ~
data$Species <- factor(data$Species, ordered = T,
                       levels = c("House",
                                  "Italian",
                                  "Spanish",
                                  NA))


# Defines the shapes to be used for each Group ~
Shapes <- as.vector(c(1, 2, 3, 13, 21, 11, 23))


# Dimensions ~ 		
# 1: D1_18.7205865386468
# 2: D2_3.05031850508455
# 3: D3_2.12431827345288


# Creates legend plot ~
MyLegend_Plot <-
  ggplot(data = data, aes_string(x = "D1_18.7205865386468", y = "D2_3.05031850508455")) +
  geom_star(aes(starshape = Population, fill = Species), size = 2.8, starstroke = .15, alpha = .85) +
  scale_fill_manual(values = c("#1E90FF", "#FFD700", "#ee0000"), na.translate = FALSE) +
  scale_starshape_manual(values = Shapes, na.translate = FALSE) +
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


# Defines the shapes to be used for each Group ~
Shapes_2 <- as.vector(c(1, 2, 3, 13, 21, 11, 23, 14))


# Combines all populations from the Faroe Islands ~
data$Species <- as.character(data$Species)
data$Population <- as.character(data$Population)
data <- data %>%
  mutate_at(c("Population", "Species"), ~replace_na(., "Target"))


# Reorders Population ~
data$Population <- factor(data$Population, ordered = T,
                          levels = c("Utrecht",
                                     "Sales",
                                     "Garderen",
                                     "Crotone",
                                     "Guglionesi",
                                     "Lesina",
                                     "Chokpak",
                                     "Target"))


# Reorders Population ~
data$Species <- factor(data$Species, ordered = T,
                       levels = c("House",
                                  "Italian",
                                  "Spanish",
                                  "Target"))


# Expands PCA_Annot by adding Labels ~
data$Labels <- ifelse(data$Species %in% c("Target"), "Target\nIndividual", "")


PCA_12_Blog <-
  ggplot(data = data, aes_string(x = "D1_18.7205865386468", y = "D2_3.05031850508455")) +
  geom_star(aes(starshape = Population, fill = Species), alpha = .85, size = 2.8, starstroke = .15) +
  scale_fill_manual(values = c("#1E90FF", "#FFD700", "#ee0000", "#d9d9d9")) +
  scale_starshape_manual(values = Shapes_2) +
  geom_label_repel(data = data, aes(label = Labels),
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
  scale_x_continuous("Dimension 1 (18.7%)",
                     #breaks = c(0.99, 1, 1.01),
                     #labels = c("0.99", "1", "1.01"),
                     #limits = c(-5, 1),
                     expand = c(.02, .02)) +
  scale_y_continuous("Dimension 2 (3.05%)",
                     #breaks = c(-.08, -.04, 0.00), 
                     #labels = c("-0.08", "-0.04", "0.00"),
                     #limits = c(-.28, .25),
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


# Isolates legend ~
MyLegendBlog <- get_legend(MyLegend_Plot)


# Gets final plot ~
PCA_Plot_Blog <- ggarrange(PCA_12_Blog, nrow = 1, legend.grob = MyLegendBlog)


# Saves plot ~
ggsave(PCA_Plot_Blog, file = "YYY.pdf",
       device = cairo_pdf, limitsize = FALSE, scale = 1, width = 10, height = 10, dpi = 600)
ggsave(PCA_Plot_Blog, file = "NLSparrow_Blog_Abridged-Ellipses--PCA.png",
       limitsize = FALSE, scale = 1, width = 9, height = 9, dpi = 600)


# Creates MDS plots ~
MDS_12 <-
ggplot(data, aes_string(x = "D1_18.7205865386468", y = "D2_3.05031850508455")) +
  geom_star(aes(starshape = Population, fill = Species), size = 2.8, alpha = .9, starstroke = .15) +
  #scale_fill_manual(values = c("#44AA99", "#F0E442", "#E69F00", "#56B4E9")) +
  scale_starshape_manual(values = Shapes) +
  #geom_mark_ellipse(aes(color = BioStatus, group = Population_cbind, filter = Population_cbind == "FaroeIslands", label = "Faroe Islands"),
  #                  label.buffer = unit(8, 'mm'), con.colour = "black", con.type = "elbow", label.fill = NA, show.legend = FALSE) +
  #geom_mark_ellipse(aes(color = BioStatus, group = Population, filter = Population == "PigeonIsland", label = "Pigeon Island"),
  #                  label.buffer = unit(40, 'mm'), con.colour = "black", con.type = "elbow", label.fill = NA, show.legend = FALSE) +
  #geom_mark_ellipse(aes(color = BioStatus, group = Population, filter = Population == "Trincomalee", label = "Trincomalee"),
  #                  label.buffer = unit(30, 'mm'), con.colour = "black", con.type = "elbow", label.fill = NA, show.legend = FALSE) +
  scale_x_continuous("Dimension 1 (18.7%)",
                     #breaks = c(-0.075, -0.05, -0.025, 0, 0.025),
                     #labels = c("-0.075", "-0.05", "-0.025", "0", "0.025"),
                     expand = c(0,0)) +
                     #limits = c(-0.073, 0.03)) +
  scale_y_continuous("Dimension 2 (3.1%)",
                     #breaks = c(-0.05, -0.025, 0, 0.025, 0.05), 
                     expand = c(0,0)) +
                     #labels = c("-0.05", "-0.025", "0", "0.025", "0.05")) +
                     #limits = c(-0.0525, 0.0525)) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.position = "right",
        legend.title = element_text(color = "#000000", size = 13),
        legend.text = element_text(size = 11),
        axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 20, face = "bold", margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text = element_text(color = "#000000", size = 13, face = "bold"),
        axis.ticks = element_line(color = "#000000", size = 0.3),
        axis.line = element_line(colour = "#000000", size = 0.3)) +
  guides(fill = guide_legend(title = "Biological Status", title.theme = element_text(size = 16, face = "bold"),
                             label.theme = element_text(size = 14),
                             override.aes = list(starshape = 21, size = 5, alpha = .9, starstroke = .15), order = 1),
         starshape = guide_legend(title = "Groups", title.theme = element_text(size = 16, face = "bold"),
                              label.theme = element_text(size = 14),
                              override.aes = list(starshape = Shapes, size = 5, alpha = .9, starstroke = .15), order = 2),
         colour = "none")


# Creates & Saves the final MDS Panel ~
ggsave(MDS_12, file = "FPG--MDS_12.pdf",
       device = cairo_pdf, scale = 1.5, width = 12, height = 8, dpi = 600)
#ggsave(MDS_12, file = "FPG--MDS_12.png",
#       scale = 1.5, width = 12, height = 8, dpi = 600)


MDS_13 <-
  ggplot(data, aes_string(x = "D1_3.18814763506453", y = "D3_1.47369056639638")) +
  geom_star(aes(starshape = Groups, fill = BioStatus), size = 2.8, alpha = .9, starstroke = .15) +
  scale_fill_manual(values = c("#44AA99", "#F0E442", "#E69F00", "#56B4E9"), labels = gsub("_", " ", levels(data$BioStatus))) +
  scale_starshape_manual(values = Shapes) +
  #geom_mark_ellipse(aes(color = BioStatus, group = Population_cbind, filter = Population_cbind == "FaroeIslands", label = "Faroe Islands"),
  #                  label.buffer = unit(8, 'mm'), con.colour = "black", con.type = "elbow", label.fill = NA, show.legend = FALSE) +
  #geom_mark_ellipse(aes(color = BioStatus, group = Population, filter = Population == "PigeonIsland", label = "Pigeon Island"),
  #                  label.buffer = unit(22, 'mm'), con.cap = unit(5, "mm"), con.colour = "black", con.type = "elbow", label.fill = NA, show.legend = FALSE) +
  #geom_mark_ellipse(aes(color = BioStatus, group = Population, filter = Population == "Trincomalee", label = "Trincomalee"),
  #                  label.buffer = unit(14, 'mm'), con.cap = unit(8.5, "mm"), con.colour = "black", con.type = "elbow", label.fill = NA, show.legend = FALSE) +
  scale_x_continuous("Dimension 1 (3.19%)",
                     breaks = c(-0.05, -0.025, 0, 0.025),
                     labels = c("-0.05", "-0.025", "0", "0.025"),
                     expand = c(0,0),
                     limits = c(-0.073, 0.03)) +
  scale_y_continuous("Dimension 3 (1.47%)",
                     breaks = c(-0.05, -0.025, 0, 0.025, 0.05), 
                     expand = c(0,0),
                     labels = c("-0.05", "-0.025", "0", "0.025", "0.05"), 
                     limits = c(-0.0525, 0.0525)) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.position = "right",
        legend.title = element_text(color = "#000000", size = 13),
        legend.text = element_text(size = 11),
        axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 20, face = "bold", margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text = element_text(color = "#000000", size = 13, face = "bold"),
        axis.ticks = element_line(color = "#000000", size = 0.3),
        axis.line = element_line(colour = "#000000", size = 0.3)) +
  guides(fill = guide_legend(title = "Biological Status", title.theme = element_text(size = 16, face = "bold"),
                             label.theme = element_text(size = 14),
                             override.aes = list(starshape = 21, size = 5, alpha = .9, starstroke = .15), order = 1),
         starshape = guide_legend(title = "Groups", title.theme = element_text(size = 16, face = "bold"),
                                  label.theme = element_text(size = 14),
                                  override.aes = list(starshape = Shapes, size = 5, alpha = .9, starstroke = .15), order = 2),
         colour = "none")


# Creates & Saves the final MDS Panel ~
ggsave(MDS_13, file = "FPG--MDS_13.pdf",
       device = cairo_pdf, scale = 1.5, width = 12, height = 8, dpi = 600)
#ggsave(MDS_13, file = "FPG--MDS_13.jpg",
#       scale = 1.5, width = 12, height = 8, dpi = 600)


MDS_23 <-
  ggplot(data, aes_string(x = "D2_1.93889781570542", y = "D3_1.47369056639638")) +
  geom_star(aes(starshape = Groups, fill = BioStatus), size = 2.8, alpha = .9, starstroke = .15) +
  scale_fill_manual(values = c("#44AA99", "#F0E442", "#E69F00", "#56B4E9"), labels = gsub("_", " ", levels(data$BioStatus))) +
  scale_starshape_manual(values = Shapes) +
  #geom_mark_ellipse(aes(color = BioStatus, group = Population_cbind, filter = Population_cbind == "FaroeIslands", label = "Faroe Islands"),
  #                  label.buffer = unit(8, 'mm'), con.colour = "black", con.type = "elbow", label.fill = NA, show.legend = FALSE) +
  #geom_mark_ellipse(aes(color = BioStatus, group = Population, filter = Population == "PigeonIsland", label = "Pigeon Island"),
  #                  label.buffer = unit(20, 'mm'), con.cap = unit(8.5, "mm"), con.colour = "black", con.type = "elbow", label.fill = NA, show.legend = FALSE) +
  #geom_mark_ellipse(aes(color = BioStatus, group = Population, filter = Population == "Trincomalee", label = "Trincomalee"),
  #                  label.buffer = unit(16, 'mm'), con.cap = unit(11.5, "mm"), con.colour = "black", con.type = "elbow", label.fill = NA, show.legend = FALSE) +
  scale_x_continuous("Dimension 2 (1.94%)",
                     breaks = c(-0.05, -0.025, 0, 0.025, 0.05),
                     labels = c("-0.05", "-0.025", "0", "0.025", "0.05"),
                     expand = c(0,0),
                     limits = c(-0.045, 0.045)) +
  scale_y_continuous("Dimension 3 (1.47%)",
                     breaks = c(-0.050, -0.025, 0, 0.025, 0.050), 
                     expand = c(0,0),
                     labels = c("-0.05", "-0.025", "0", "0.025", "0.05"), 
                     limits = c(-0.0525, 0.0525)) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.position = "right",
        legend.title = element_text(color = "#000000", size = 13),
        legend.text = element_text(size = 11),
        axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 20, face = "bold", margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text = element_text(color = "#000000", size = 13, face = "bold"),
        axis.ticks = element_line(color = "#000000", size = 0.3),
        axis.line = element_line(colour = "#000000", size = 0.3)) +
  guides(fill = guide_legend(title = "Biological Status", title.theme = element_text(size = 16, face = "bold"),
                             label.theme = element_text(size = 14),
                             override.aes = list(starshape = 21, size = 5, alpha = .9, starstroke = .15), order = 1),
         starshape = guide_legend(title = "Groups", title.theme = element_text(size = 16, face = "bold"),
                                  label.theme = element_text(size = 14),
                                  override.aes = list(starshape = Shapes, size = 5, alpha = .9, starstroke = .15), order = 2),
         colour = "none")


# Creates & Saves the final MDS Panel ~
ggsave(MDS_23, file = "FPG--MDS_23.pdf",
       device = cairo_pdf, width = 12, height = 8, scale = 1.5, dpi = 600)
#ggsave(MDS_23, file = "FPG--MDS_23.jpg",
#       width = 12, height = 8, scale = 1.5, dpi = 600)


# Creates MDS_SI ~
MDS_SI <- ggarrange(MDS_13, MDS_23, labels = c("A", "B"), ncol = 1, nrow = 2,
          common.legend = TRUE, legend = "right")


# Saves MDS_SI ~
ggsave(MDS_SI, file = "FPG--MDS_SI.pdf",
       device = cairo_pdf, width = 12, height = 14, scale = 1.5, dpi = 600)
#ggsave(MDS_SI, file = "FPG--MDS_SI.jpg",
#       width = 12, height = 14, scale = 1.5, dpi = 600)


#
##
### The END ~~~~~