### The BEGINNING ~~~~~
##
# ~ Creates NLSparrow--GlobalDepth | By George Pacheco.


# Cleans the environment ~ 
rm(list=ls())


# Sets working directory ~
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Loads required packages ~
pacman::p_load(scales, extrafont, tidyverse, reshape2, lemon)


# Imports extra fonts ~
loadfonts(device = "win", quiet = TRUE)


# Loads coverage data ~
lgz <- list()
m <- list()
listgz <- dir(pattern = ".gz")
for (k in 1:length(listgz)){
  lgz[[k]] <- read.table(gzfile(listgz[k]))
  colnames(lgz[[k]]) <- c("CHR", "Start", "End", "Depth")
  lgz[[k]]$Species <- gsub("*NLSparrow_MinInd90.GlobalDepth.Subset.gz", "", listgz[k])
  lgz[[k]]$Species <- gsub("BSG_", "", lgz[[k]]$Species)}


# Melts PSMC data ~
fulldfUp <- melt(lgz, id = c("CHR", "Start", "End", "Depth", "DataType"))
fulldfUp$Type <-  ""


# Corrects Population names ~
levels(fulldfUp$Species <- sub("EuropeanFlounder", "European Flounder", fulldfUp$Species))
levels(fulldfUp$Species <- sub("AtlanticCod", "Atlantic Cod", fulldfUp$Species))
levels(fulldfUp$Species <- sub("AtlanticHerring", "Atlantic Herring", fulldfUp$Species))


# Reorders Population ~
fulldfUp$Species <- factor(fulldfUp$Species, ordered = T,
                              levels = c("Turbot",
                                         "European Flounder",
                                         "Atlantic Cod",
                                         "Atlantic Herring",
                                         "Lumpfish"))


# Subsamples data ~
Turbot <- subset(fulldfUp, Species == "Turbot")
EuropeanFlounder <- subset(fulldfUp, Species == "European Flounder")
AtlanticCod <- subset(fulldfUp, Species == "Atlantic Cod")
AtlanticHerring <- subset(fulldfUp, Species == "Atlantic Herring")
Lumpfish <- subset(fulldfUp, Species == "Lumpfish")

quantile(fulldfUp$Depth, c(.975))


# Creates the plot ~
GlobalCoverage <-
 ggplot(fulldfUp, aes(x = Depth, fill = Type, colour = Type)) +
  geom_density(alpha = .15, size = .3) +
  geom_vline(data = fulldfUp, aes(xintercept = quantile(Depth, c(.975))), colour = "#df65b0") +
  facet_rep_grid(Species ~. , scales = "free_x") +
  scale_fill_manual(values = c("#df65b0")) +
  scale_colour_manual(values = c("#df65b0")) +
  scale_x_continuous("Global Depth (X)",
                     breaks = c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200),
                     labels = c("100", "200", "300", "400", "500", "600", "700", "800", "900", "1000", "1100", "1200"),
                     limits = c(0, 1215),
                     expand = c(0, 0)) +
  scale_y_continuous("Density",
                     breaks = c(.001, .002, .003, .004, .005),
                     labels = c("0.001", "0.002", "0.003", "0.004", "0.005"),
                     expand = c(0, 0)) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.grid.major = element_line(color = "#ededed", linetype = "dashed", linewidth = .00005),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.spacing.y = unit(1, "cm"),
        axis.title.x = element_text(size = 15, face = "bold", color = "#000000", margin = margin(t = 25, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 15, face = "bold", color = "#000000", margin = margin(t = 0, r = 25, b = 0, l = 0)),
        axis.ticks = element_line(linewidth = .3, color = "#000000"),
        strip.text = element_text(colour = "#000000", size = 11, face = "bold"),
        strip.background.y = element_rect(colour = "#000000", fill = "#d6d6d6", linewidth = .3),
        axis.text.x = element_text(size = 9, color = "#000000", face = "bold", angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 9, color = "#000000", face = "bold"),
        axis.line = element_line(colour = "#000000", linewidth = .3),
        legend.position = "none")

  
# Saves plot ~
ggsave(GlobalCoverage, file = "BSG_Combined--GlobalDepth.pdf",
       scale = 1, width = 12, height = 12, device = cairo_pdf, dpi = 600)
ggsave(GlobalCoverage, file = "BSG_Combined--GlobalDepth.jpg",
       scale = 1, width = 12, height = 12, dpi = 600)


#
##
### The END ~~~~~