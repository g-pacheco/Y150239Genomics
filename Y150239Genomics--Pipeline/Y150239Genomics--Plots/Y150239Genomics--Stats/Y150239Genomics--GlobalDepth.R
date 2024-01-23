### The BEGINNING ~~~~~
##
# ~ Creates Y150239Genomics--GlobalDepth | By George Pacheco.


# Cleans the environment ~ 
rm(list=ls())


# Sets working directory ~
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Loads required packages ~
pacman::p_load(scales, extrafont, tidyverse, reshape2, lemon)


# Imports extra fonts ~
font_import()


# Loads data ~
fulldf <- read.table("AllSamples_haplotypecaller.raw.vcf.Filtered.NoKinship.NoTree.MeanDepth.ldepth.mean", header = TRUE)
fulldf$Type <- ""


# Expands fulldf by adding chrtype ~
fulldf$chrtype <- ifelse(grepl("chrZ", fulldf$CHROM), "Allosome (Z)",  "Autosomes")


# Reorders chrtype ~
fulldf$chrtype <- factor(fulldf$chrtype, ordered = T,
                         levels = c("Autosomes",
                                    "Allosome (Z)"))


# Get quantile values ~
quantiles_df <- fulldf %>%
  group_by(chrtype) %>%
  summarize(x_quantile = quantile(MEAN_DEPTH, .95))


# Expands quantiles_df by adding annotation ~
quantiles_df <- quantiles_df %>%
  add_column(group = c("Autosomes", "Allosome (Z)")) %>%
  add_column(label = c(sprintf("Quantile 95%%: %.0fX", quantiles_df$x_quantile))) %>%
  add_column(value_x = c(62.5, 62.5)) %>%
  add_column(value_y = c(1.25, .035)) %>%
  add_column(Type = "")


# Custom y-axis breaks ~
breaks_fun <- function(y){
  caseVal <- max(y)
  if (caseVal > 1){
    c(.5, 1, 1.5)}
  else {
    c(.01, .02, .03, .04)}}


# Custom y-axis limits ~
limits_fun <- function(x){
  limitVal <- max(x)
  if (limitVal > 1){
    c(0, 1.575)}
  else {
    c(0, .0425)}}


# Creates the plot ~
GlobalCoverage <-
 ggplot(fulldf, aes(x = MEAN_DEPTH, fill = Type, colour = Type), alpha = .1) +
  geom_density(alpha = .15, linewidth = .3, adjust = 1) +
  facet_rep_grid(chrtype ~., scales = "free") +
  scale_fill_manual(values = c("#fa9fb5")) +
  scale_colour_manual(values = c("#000000")) +
  geom_label(data = quantiles_df, aes(x = value_x, y = value_y, label = label), alpha = .1,
             size = 3.7, fontface = "bold", family = "Optima", show.legend = FALSE) +
  scale_x_continuous("Mean Depth (X)",
                     breaks = c(25, 50, 75, 100),
                     labels = c("25", "50", "75", "100"),
                     limits = c(0, 102),
                     expand = c(0, 0)) +
  scale_y_continuous("Density",
                     breaks = breaks_fun,
                     limits = limits_fun,
                     expand = c(0, 0)) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.grid.major = element_line(color = "#d9d9d9", linetype = "dashed", linewidth = .05),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.spacing.y = unit(1, "cm"),
        axis.title.x = element_text(size = 16, face = "bold", color = "#000000", margin = margin(t = 25, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 16, face = "bold", color = "#000000", margin = margin(t = 0, r = 25, b = 0, l = 0)),
        axis.ticks = element_line(linewidth = .3, color = "#000000"),
        axis.text.x = element_text(size = 10, color = "#000000", face = "bold", angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 10, color = "#000000", face = "bold"),
        strip.text = element_text(colour = "#000000", size = 13, face = "bold", family = "Optima"),
        axis.line = element_line(colour = "#000000", linewidth = .3),
        strip.background.y = element_rect(colour = "#000000", fill = "#d6d6d6", linewidth = .3),
        legend.position = "none")

  
# Saves plot ~
ggsave(GlobalCoverage, file = "Y150239Genomics--MeanDepth.pdf",
       scale = 1, width = 9, height = 9, device = cairo_pdf, dpi = 600)
ggsave(GlobalCoverage, file = "Y150239Genomics--MeanDepth.jpg",
       scale = 1, width = 9, height = 9, dpi = 600)


#
##
### The END ~~~~~