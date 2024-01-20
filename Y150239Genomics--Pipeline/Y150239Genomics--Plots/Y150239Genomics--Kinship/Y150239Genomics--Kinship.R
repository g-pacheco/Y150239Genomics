### The BEGINNING ~~~~~
##
# Y150239Genomics--Kinship by George Pacheco ~


# Cleans the environment ~ 
rm(list=ls())


# Sets working directory ~
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Loads required packages ~
pacman::p_load(optparse, tidyverse, plyr, RColorBrewer, extrafont, ggforce, ggstar, ggrepel, RcppCNPy, reshape2, lemon, plotly,
               gridExtra, grid, cowplot, patchwork, ggpubr, rphylopic, viridis)


# Imports data while incorporating annotation ~
annot <- list()
rab <- list()
annotL <- dir(pattern = ".labels")
rabL <- dir(pattern = ".res")
for (k in 1:length(annotL)){
  annot[[k]] <- read.table(annotL[k], sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  colnames(annot[[k]]) <- c("Ind1")
  annot[[k]]$Ind2 <- c(annot[[k]]$Ind1[-c(1)], NA)
  rab[[k]] <- read.table(rabL[k], header = TRUE, stringsAsFactors = FALSE)
  rab[[k]]$a <- rab[[k]]$a + 1
  rab[[k]]$Population <- gsub("AllSamples_haplotypecaller.raw.vcf.Filtered.", "", rabL[k])
  rab[[k]]$Population <- gsub(".res", "", rab[[k]]$Population)
  rab[[k]]$a <- annot[[k]]$Ind1[match(rab[[k]]$a, seq_along(annot[[k]]$Ind1))]
  rab[[k]]$b <- annot[[k]]$Ind2[match(rab[[k]]$b, seq_along(annot[[k]]$Ind2))]
  rab[[k]]$Pair <- paste(rab[[k]]$a,"Vs",rab[[k]]$b)
  rab[[k]] <- rab[[k]] %>% select(a, b, rab, Pair, Population)}


# Expands list ~
fulldf <- bind_rows(rab)


# Reorders Population ~
fulldf$Population <- factor(fulldf$Population, ordered = T,
                            levels = c("Utrecht",
                                       "Sales",
                                       "Crotone",
                                       "Guglionesi", 
                                       "Lesina",
                                       "Chokpak"))


# Creates plot (Boxplot) ~
Kinship_Plot_Boxplot <-
  ggplot(fulldf, aes(x = Population, y = rab)) +
  geom_boxplot(fill = "#ffffff", colour = "#000000", show.legend = FALSE, alpha = .9, linewidth = .3, width = .15) +
  geom_point(data = subset(fulldf, rab >= .125), color = "#df65b0") +
  geom_hline(yintercept = .125, linetype = "twodash", color = "#df65b0", linewidth = .3) +
  geom_label_repel(data = subset(fulldf, rab >= .125), aes(label = Pair),
                   size = 3, nudge_x = .1, fontface = "bold", family = "Times New Roman") +
  scale_x_discrete(expand = c(.05, .05)) +
  scale_y_continuous ("Rab", 
                      breaks = c(.1, .2, .3, .4, .5), 
                      labels = c("0.1", "0.2", "0.3", "0.4", "0.5"),
                      limits = c(0, .5025),
                      expand = c(.005, .005)) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.margin = margin(t = 0, b = 0, r = 15, l = 15),
        legend.box = "vertical",
        legend.box.margin = margin(t = 20, b = 30, r = 0, l = 0),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "#000000", size = 12, face = "bold", margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.x = element_text(color = "#000000", size = 10, face = "bold", angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(color = "#000000", size = 9, face = "bold"),
        axis.ticks = element_line(color = "#000000", linewidth = .3),
        strip.text = element_text(colour = "#000000", size = 13, face = "bold"),
        strip.background = element_rect(colour = "#000000", fill = "#d6d6d6", linewidth = .3),
        axis.line = element_line(colour = "#000000", linewidth = .3)) +
  guides(fill = guide_legend(title = "Rab", title.theme = element_text(size = 16, face = "bold"),
                             label.theme = element_text(size = 15)))


# Saves plot (Boxplot) ~
ggsave(Kinship_Plot_Boxplot, file = "Y150239Genomics--Kinship_Boxplot.pdf",
       device = cairo_pdf, limitsize = FALSE, scale = 1, width = 8, height = 5, dpi = 600)
ggsave(Kinship_Plot_Boxplot, file = "Y150239Genomics--Kinship_Boxplot.jpeg",
      limitsize = FALSE, scale = 1, width = 8, height = 5, dpi = 600)


# Creates plot (Heatmap) ~
Kinship_Plot_Heatmap <-
  ggplot(fulldf, aes(a, b, fill = rab)) + 
  geom_tile(colour = "#ffffff") +
  scale_fill_continuous(trans = "reverse") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(limits = rev, expand = c(0, 0)) +
  facet_wrap(Population ~., scales = "free", ncol = 2) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
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
        strip.text = element_text(colour = "#000000", size = 13, face = "bold"),
        strip.background = element_rect(colour = "#000000", fill = "#d6d6d6", linewidth = .3),
        axis.line = element_line(colour = "#000000", linewidth = .3)) +
  guides(fill = guide_legend(title = "Rab", title.theme = element_text(size = 16, face = "bold"),
                             label.theme = element_text(size = 15)))


# Saves plot (Heatmap) ~
ggsave(Kinship_Plot_Heatmap, file = "Y150239Genomics--Kinship_Heatmap.pdf",
       device = cairo_pdf, limitsize = FALSE, scale = 1, width = 12, height = 12, dpi = 600)
ggsave(Kinship_Plot_Heatmap, file = "Y150239Genomics--Kinship_Heatmap.jpeg",
       limitsize = FALSE, scale = 1, width = 8, height = 5, dpi = 600)


#
##
### The END ~~~~~


# Melts data ~
pops <- union(Kinship$a, Kinship$b)
n <- length(pops)


# Creates Kinship-Kinship matrix ~
KinshipUp <- matrix(0, nrow = n, ncol = n, dimnames = list(pops, pops))
for (i in 1:nrow(Kinship)) {
  KinshipUp[Kinship[i, "a"], Kinship[i, "b"]] = Kinship[i, "rab"]
  KinshipUp[Kinship[i, "b"], Kinship[i, "a"]] = Kinship[i, "rab"]}


# Melts matrix ~
fulldf <- melt(KinshipUp)