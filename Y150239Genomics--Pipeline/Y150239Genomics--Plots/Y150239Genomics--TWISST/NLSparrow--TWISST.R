### The BEGINNING ~~~~~
##
# NLSparrow--TWISST | by George Pacheco ~


# Cleans the environment ~ 
rm(list=ls())


# Sets working directory ~
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Loads packages ~
pacman::p_load(tidyverse, scales, reshape2, lemon)


# Imports weights ~
Wlgz <- list()
#Wlistgz <- dir(path = "/../../../../LargeFiles/", pattern = ".csv.gz")
Wlistgz <- dir(pattern = ".gz")
for (k in 1:length(Wlistgz)){
  Wlgz[[k]] <- read.table(gzfile(Wlistgz[k]))[-(1:1), ]
  colnames(Wlgz[[k]]) <- c("Target ~ Spanish", "Target ~ Italian", "Target ~ House")
  Wlgz[[k]]$CHR <- gsub("AllSamples_haplotypecaller.raw.vcf.Filtered.MAF20.", "", Wlistgz[k])
  Wlgz[[k]]$CHR <- gsub(".SW250.Weights.csv.gz", "", Wlgz[[k]]$CHR)}


# Melts weights ~
WeightsDF <- reshape2::melt(Wlgz)
WeightsDF <- WeightsDF[, -ncol(WeightsDF)]


# Imports windows´ data ~
Wilgz <- list()
#Wilistgz <- dir(path = "../../../../LargeFiles/", pattern = ".txt")
Wilistgz <- dir(pattern = ".txt")
for (k in 1:length(Wilistgz)){
  Wilgz[[k]] <- read.table(Wilistgz[k])[-1, ]
  colnames(Wilgz[[k]]) <- c("Scaffold", "Start", "End", "Mid", "Sites", "lnL")}


# Melts windows ~
WindowsDF <- reshape2::melt(Wilgz)
WindowsDF <- WindowsDF[, -ncol(WindowsDF)]


# Merges WeightsDF & WindowsDF ~
fulldf <- cbind(WeightsDF, WindowsDF)


# Converts DF from wide into long ~
fulldfUp <- gather(fulldf, Weight, Value, "Target ~ Spanish", "Target ~ Italian", "Target ~ House")


# Reorders CHR ~
fulldfUp$CHR <- factor(fulldfUp$CHR, ordered = T,
                       levels = c("chr1", "chr1A", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr17",
                                  "chr18", "chr19", "chr20", "chr21", "chr22", "chr23", "chr24", "chr26", "chr27", "chr28"))


# Reorders Weight ~
fulldfUp$Weight <- factor(fulldfUp$Weight, ordered = T,
                       levels = c("Target ~ House",
                                  "Target ~ Italian",
                                  "Target ~ Spanish"))


# Corrects the y-strip facet labels ~
y_strip_labels <- c("chr1" = "CHR 01", "chr1A" = "CHR 01A", "chr2" = "CHR 02", "chr3" = "CHR 03", "chr4" = "CHR 04", "chr5" = "CHR 05", "chr6" = "CHR 06", "chr7" = "CHR 07",
                    "chr8" = "CHR 08", "chr9" = "CHR 09", "chr10" = "CHR 10", "chr11" = "CHR 11", "chr12" = "CHR 12", "chr13" = "CHR 13", "chr14" = "CHR 14", "chr15" = "CHR 15",
                    "chr17" = "CHR 17", "chr18" = "CHR 18", "chr19" = "CHR 19", "chr20" = "CHR 20", "chr21" = "CHR 21", "chr22" = "CHR 22", "chr23" = "CHR 23",
                    "chr24" = "CHR 24", "chr26" = "CHR 26", "chr27" = "CHR 27", "chr28" = "CHR 28")


# Creates plot ~
Weights_Area_Real_Plot <-
  ggplot() +
  geom_area(data = fulldfUp, aes(x = as.numeric(Mid), y = as.numeric(Value), fill = Weight, group = Weight),
            position = "identity", colour = "#000000", alpha = .3, linetype = 1, linewidth = .2) +
  facet_rep_grid(CHR ~., scales = "free_y", labeller = labeller(CHR = y_strip_labels)) +
  scale_x_continuous("Genomic Position",
                     breaks = c(25000000, 50000000, 75000000, 100000000, 125000000), 
                     labels = c("25Mb", "50Mb", "75Mb", "100Mb", "125Mb"),
                     limits = c(0, 147697000),
                     expand = c(0, 0)) +
  scale_y_continuous("Weights Across Chrmosomes",
                     breaks = c(50000, 100000, 150000, 200000), 
                     labels = c("50K", "100K", "150K", "200K"),
                     limits = c(0, 201000),
                     expand = c(0, 0)) +
  scale_fill_manual(values = c("#c2a5cf", "#3288bd", "#d53e4f")) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.border = element_blank(),
        panel.grid.major = element_line(color = "#ededed", linetype = "dashed", linewidth = .00005),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "#000000", linewidth = .3),
        axis.title.x = element_text(size = 20, face = "bold", color = "#000000", margin = margin(t = 30, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 20, face = "bold", color = "#000000", margin = margin(t = 0, r = 30, b = 0, l = 0)),
        axis.text = element_text(colour = "#000000", size = 8, face = "bold"),
        axis.ticks = element_line(color = "#000000", linewidth = .3),
        strip.background.y = element_rect(colour = "#000000", fill = "#d6d6d6", linewidth = .3),
        strip.text = element_text(colour = "#000000", size = 11.5, face = "bold"),
        legend.position = "top",
        legend.margin = margin(t = 0, b = 0, r = 0, l = 0),
        legend.box.margin = margin(t = 30, b = 40, r = 0, l = 0),
        legend.key = element_rect(fill = NA),
        legend.background = element_blank()) +
  guides(colour = guide_legend(title = "Topologies", title.theme = element_text(size = 21, face = "bold"),
                               label.theme = element_text(size = 19), override.aes = list(linewidth = 1.75, linetype = 1)),
         fill = guide_legend(title = "Topologies", title.theme = element_text(size = 21, face = "bold"),
                             label.theme = element_text(size = 19), override.aes = list(linewidth = .3, linetype = 1)))


# Saves plot ~
ggsave(Weights_Area_Real_Plot, file = "NLSparrow--TWISST_AreaReal.pdf",
       device = cairo_pdf, limitsize = FALSE, width = 30, height = 30, scale = 1, dpi = 100)
ggsave(Weights_Area_Real_Plot, file = "NLSparrow--TWISST_AreaReal.jpeg",
       limitsize = FALSE, width = 30, height = 35, scale = 1, dpi = 100)


# Creates plot ~
Weights_Area_Real_Smooth_Plot <-
  fulldfUp %>%
  ggplot(aes(x = as.numeric(Mid), y = as.numeric(Value), fill = Weight, group = Weight)) +
  stat_smooth(geom = 'area', method = 'loess', position = "identity", linetype = 1, linewidth = .2, colour = "#000000", span = .15, alpha = .3) +
  facet_rep_grid(CHR ~., scales = "free_y", labeller = labeller(CHR = y_strip_labels)) +
  scale_x_continuous("Genomic Position",
                     breaks = c(25000000, 50000000, 75000000, 100000000, 125000000), 
                     labels = c("25Mb", "50Mb", "75Mb", "100Mb", "125Mb"),
                     limits = c(0, 147697000),
                     expand = c(0, 0)) +
  scale_y_continuous("Weights Across Chrmosomes",
                     breaks = c(50000, 100000, 150000, 200000), 
                     labels = c("50K", "100K", "150K", "200K"),
                     limits = c(0, 201000),
                     expand = c(0, 0)) +
  scale_fill_manual(values = c("#c2a5cf", "#3288bd", "#d53e4f")) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.border = element_blank(),
        panel.grid.major = element_line(color = "#ededed", linetype = "dashed", linewidth = .00005),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "#000000", linewidth = .3),
        axis.title.x = element_text(size = 20, face = "bold", color = "#000000", margin = margin(t = 30, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 20, face = "bold", color = "#000000", margin = margin(t = 0, r = 30, b = 0, l = 0)),
        axis.text = element_text(colour = "#000000", size = 8, face = "bold"),
        axis.ticks = element_line(color = "#000000", linewidth = .3),
        strip.background.y = element_rect(colour = "#000000", fill = "#d6d6d6", linewidth = .3),
        strip.text = element_text(colour = "#000000", size = 11.5, face = "bold"),
        legend.position = "top",
        legend.margin = margin(t = 0, b = 0, r = 0, l = 0),
        legend.box.margin = margin(t = 30, b = 40, r = 0, l = 0),
        legend.key = element_rect(fill = NA),
        legend.background = element_blank()) +
  guides(colour = guide_legend(title = "Topologies", title.theme = element_text(size = 21, face = "bold"),
                               label.theme = element_text(size = 19), override.aes = list(linewidth = 1.75, linetype = 1)),
         fill = guide_legend(title = "Topologies", title.theme = element_text(size = 21, face = "bold"),
                             label.theme = element_text(size = 19), override.aes = list(linewidth = .3, linetype = 1)))


# Saves plot ~
ggsave(Weights_Area_Real_Smooth_Plot, file = "NLSparrow--TWISST_AreaRealSmooth.pdf",
       device = cairo_pdf, limitsize = FALSE, width = 30, height = 30, scale = 1, dpi = 100)
ggsave(Weights_Area_Real_Smooth_Plot, file = "NLSparrow--TWISST_AreaRealSmooth.jpeg",
       limitsize = FALSE, width = 30, height = 35, scale = 1, dpi = 100)

  
# Creates plot ~
Weights_Area_Fill_Plot <-
  ggplot() +
  geom_area(data = fulldfUp, aes(x = as.numeric(Mid), y = as.numeric(Value), fill = Weight, group = Weight),
            position = "fill", colour = "#000000", alpha = .3, linetype = 1, linewidth = .2) +
  facet_rep_grid(CHR ~., scales = "free_y", labeller = labeller(CHR = y_strip_labels)) +
  scale_x_continuous("Genomic Position",
                   breaks = c(25000000, 50000000, 75000000, 100000000, 125000000), 
                   labels = c("25Mb", "50Mb", "75Mb", "100Mb", "125Mb"),
                   limits = c(0, 147697000),
                   expand = c(0, 0)) +
  scale_y_continuous("Weights Across Chrmosomes",
                    breaks = c(.25, .50, .75), 
                    labels = c("25%", "50%", "75%"),
                    expand = c(0, 0)) +
  scale_fill_manual(values = c("#c2a5cf", "#3288bd", "#d53e4f")) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.border = element_blank(),
        panel.grid.major = element_line(color = "#ededed", linetype = "dashed", linewidth = .00005),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "#000000", linewidth = .3),
        axis.title.x = element_text(size = 20, face = "bold", color = "#000000", margin = margin(t = 30, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 20, face = "bold", color = "#000000", margin = margin(t = 0, r = 30, b = 0, l = 0)),
        axis.text = element_text(colour = "#000000", size = 8, face = "bold"),
        axis.ticks = element_line(color = "#000000", linewidth = .3),
        strip.background.y = element_rect(colour = "#000000", fill = "#d6d6d6", linewidth = .3),
        strip.text = element_text(colour = "#000000", size = 11.5, face = "bold"),
        legend.position = "top",
        legend.margin = margin(t = 0, b = 0, r = 0, l = 0),
        legend.box.margin = margin(t = 30, b = 40, r = 0, l = 0),
        legend.key = element_rect(fill = NA),
        legend.background = element_blank()) +
  guides(colour = guide_legend(title = "Topologies", title.theme = element_text(size = 21, face = "bold"),
                               label.theme = element_text(size = 19), override.aes = list(linewidth = 1.75, linetype = 1)),
         fill = guide_legend(title = "Topologies", title.theme = element_text(size = 21, face = "bold"),
                             label.theme = element_text(size = 19), override.aes = list(linewidth = .3, linetype = 1)))


# Saves plot ~
ggsave(Weights_Area_Fill_Plot, file = "NLSparrow--TWISST_AreaFill.pdf",
       device = cairo_pdf, limitsize = FALSE, width = 30, height = 30, scale = 1, dpi = 100)
ggsave(Weights_Area_Fill_Plot, file = "NLSparrow--TWISST_AreaFill.jpeg",
       limitsize = FALSE, width = 30, height = 35, scale = 1, dpi = 100)


# Creates plot ~
Weights_Line_Plot <-
  ggplot() +
  geom_line(data = fulldfUp, aes(x = as.numeric(Mid), y = as.numeric(Value), colour = Weight, group = Weight),
            position = "identity", linetype = 1, linewidth = .2) +
  facet_rep_grid(CHR ~., scales = "free_y", labeller = labeller(CHR = y_strip_labels)) +
  scale_x_continuous("Genomic Position",
                     breaks = c(25000000, 50000000, 75000000, 100000000, 125000000), 
                     labels = c("25Mb", "50Mb", "75Mb", "100Mb", "125Mb"),
                     limits = c(0, 147697000),
                     expand = c(0, 0)) +
  scale_y_continuous("Weights Across Chrmosomes",
                     breaks = c(50000, 100000, 150000, 200000), 
                     labels = c("50K", "100K", "150K", "200K"),
                     limits = c(0, 201000),
                     expand = c(0, 0)) +
  scale_colour_manual(values = c("#c2a5cf", "#3288bd", "#d53e4f")) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.border = element_blank(),
        panel.grid.major = element_line(color = "#ededed", linetype = "dashed", linewidth = .00005),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "#000000", linewidth = .3),
        axis.title.x = element_text(size = 20, face = "bold", color = "#000000", margin = margin(t = 30, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 20, face = "bold", color = "#000000", margin = margin(t = 0, r = 30, b = 0, l = 0)),
        axis.text = element_text(colour = "#000000", size = 8, face = "bold"),
        axis.ticks = element_line(color = "#000000", linewidth = .3),
        strip.background.y = element_rect(colour = "#000000", fill = "#d6d6d6", linewidth = .3),
        strip.text = element_text(colour = "#000000", size = 11.5, face = "bold"),
        legend.position = "top",
        legend.margin = margin(t = 0, b = 0, r = 0, l = 0),
        legend.box.margin = margin(t = 30, b = 40, r = 0, l = 0),
        legend.key = element_rect(fill = NA),
        legend.background = element_blank()) +
  guides(colour = guide_legend(title = "Topologies", title.theme = element_text(size = 21, face = "bold"),
                               label.theme = element_text(size = 19), override.aes = list(linewidth = 1.75, linetype = 1)),
         fill = guide_legend(title = "Topologies", title.theme = element_text(size = 21, face = "bold"),
                             label.theme = element_text(size = 19), override.aes = list(linewidth = .3, linetype = 1)))


# Saves plot ~
ggsave(Weights_Line_Plot, file = "NLSparrow--TWISST_LineReal.pdf",
       device = cairo_pdf, limitsize = FALSE, width = 30, height = 30, scale = 1, dpi = 100)
ggsave(Weights_Line_Plot, file = "NLSparrow--TWISST_LineReal.jpeg",
       limitsize = FALSE, width = 30, height = 35, scale = 1, dpi = 100)


#
##
### The END ~~~~~