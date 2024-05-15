### The BEGINNING ~~~~~
##
# ~ Plots NLSparrow--evalAdmix | First written by Jose Samaniego with later modifications by George Pacheco.


# Cleans the environment ~ 
rm(list=ls())


# Sets working directory ~
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Loads required packages ~
pacman::p_load(tidyverse)
#chattr::chattr_app()
source("visFuns.R")

# Imports data while incorporating annotation ~
corres <- list()
annot <- list()
corres_files <- dir(pattern = ".corres")
annot_files <- dir(pattern = ".labels")
for (k in seq_along(annot_files)) {
  annot[[k]] <- read.table(annot_files[k], sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  colnames(annot[[k]]) <- c("Annot")
  annot[[k]]$Population_1 <- ifelse(grepl("FR0", annot[[k]]$Annot), "Sales",
                             ifelse(grepl("KAZ", annot[[k]]$Annot), "Chokpak",
                             ifelse(grepl("Lesina", annot[[k]]$Annot), "Lesina",
                             ifelse(grepl("Crotone", annot[[k]]$Annot), "Crotone",
                             ifelse(grepl("Guglionesi", annot[[k]]$Annot), "Guglionesi",
                             ifelse(grepl("PI22NLD0001M", annot[[k]]$Annot), "Y150239",
                             ifelse(grepl("PD22NLD0146F", annot[[k]]$Annot), "Garderen",
                             ifelse(grepl("PD22NLD0147F", annot[[k]]$Annot), "Garderen",
                             ifelse(grepl("PDOM2022NLD0077M", annot[[k]]$Annot), "Meerkerk",
                             ifelse(grepl("PDOM2022NLD0", annot[[k]]$Annot), "Utrecht", "Error"))))))))))
  corres_df <- as.data.frame(read.table(corres_files[k]))
  rownames(corres_df) <- annot[[k]]$Annot
  colnames(corres_df) <- annot[[k]]$Annot
  corres_df$CHRType <- str_extract(corres_files[k], "(Allosome|Autosomes)")
  corres_df$K <- str_extract(corres_files[k], "(K2|K3|K4|K5|K6|K7)")
  corres_df$K <- ifelse(grepl("K2", corres_df$K), "K = 2",
                 ifelse(grepl("K3", corres_df$K), "K = 3",
                 ifelse(grepl("K4", corres_df$K), "K = 4",
                 ifelse(grepl("K5", corres_df$K), "K = 5",
                 ifelse(grepl("K6", corres_df$K), "K = 6",
                 ifelse(grepl("K7", corres_df$K), "K = 7", "Error"))))))
  corres[[k]] <- corres_df}


# Defines compute_individual_correlations function ~
compute_individual_correlations <- function(corres_ind_list, pop_list, ord_list) {
  result <- lapply(1:length(corres_ind_list), function(i) {
    corres_ind <- corres_ind_list[[i]]
    pop <- pop_list[[i]]
    ord <- ord_list[[i]]
    annotation_df <- data.frame(Sample_ID_1 = rownames(corres_ind), Population_1 = pop)
    for (j in 1:length(pop)) {
      pop_name <- pop[j]
      if (pop_name %in% names(ord)) {
        pop_ord <- ord[[pop_name]]
        chr_type <- unique(pop_ord$CHRType)
        k_val <- unique(pop_ord$K)
        annotation_df$CHRType[pop == pop_name] <- chr_type
        annotation_df$K[pop == pop_name] <- k_val}}
    # Replace NA values with 10
    # corres_ind[is.na(corres_ind)] <- 10
    
    indiv_matrix <- cbind(annotation_df, corres_ind)
    return(indiv_matrix)})
  return(result)}


# Applies compute_on_matrices function ~
ind_list <- compute_individual_correlations(corres_ind_list = corres,
                                            pop_list = lapply(annot, function(ann) ann$Population_1),
                                            ord_list = lapply(corres, orderInds))


# Removes non-computable columns from main_list ~
exclude_cols <- as.vector(c("Sample_ID_1", "CHRType", "K", "Population_1"))
ind_list_red <- lapply(ind_list, function(cor_mat) {
                       cor_mat[, !(colnames(cor_mat) %in% exclude_cols)]})


# Defines vector of populations ~
pop <- as.vector(annot[[1]]$Population_1)


# Defines compute_mean_correlations function ~
compute_mean_correlations <- function(cor_mat_list, pop) {
  unique_pops <- unique(pop)
  num_pops <- length(unique_pops)
  
  process_single_matrix <- function(cor_mat) {
    mean_cor_df <- data.frame(matrix(ncol = num_pops, nrow = num_pops))
    rownames(mean_cor_df) <- unique_pops
    colnames(mean_cor_df) <- unique_pops
    
    # Add Population_1 column
    mean_cor_df$Population_1 <- unique(pop)
    
    for (i1 in 1:num_pops) {
      for (i2 in 1:num_pops) {
        p1 <- unique_pops[i1]
        p2 <- unique_pops[i2]
        
        # Ensure the subsetting is correct
        indices_p1 <- which(pop == p1)
        indices_p2 <- which(pop == p2)
        
        cor_values <- cor_mat[indices_p1, indices_p2]
        mean_cor_df[p1, p2] <- mean(cor_values[!is.na(cor_values)])}}
    
    # Replace NaN values with 0
    mean_cor_df[is.na(mean_cor_df)] <- 0
    
    return(mean_cor_df)}
  
  mean_cor_list <- lapply(cor_mat_list, process_single_matrix)
  return(mean_cor_list)}


# Applies compute_mean_correlation function ~
mean_list <- compute_mean_correlations(ind_list_red, pop)


# Defines annotate_mean_df function ~
annotate_mean_df <- function(larger_list, smaller_list) {
  # Iterate through each pair of data frames
  for (i in seq_along(larger_list)) {
    larger_df <- larger_list[[i]]
    smaller_df <- smaller_list[[i]]
    
    # Extract CHRType and K columns from larger data frame
    CHRType <- larger_df$CHRType[1]
    K <- larger_df$K[1]
    
    # Add CHRType and K columns to smaller data frame
    smaller_df$CHRType <- CHRType
    smaller_df$K <- K
    
    # Assign the modified smaller data frame back to the list
    smaller_list[[i]] <- smaller_df}
  return(smaller_list)}


# Applies annotate_mean_df function ~
mean_list <- annotate_mean_df(ind_list, mean_list)


# Convert all matrices to numeric type
ind_list <- lapply(ind_list, as.matrix)
ind_list <- lapply(ind_list, as.data.frame)
mean_list <- lapply(mean_list, as.matrix)
mean_list <- lapply(mean_list, as.data.frame)


# Combine all matrices into a single dataset
ind_combined <- bind_rows(ind_list, .id = "Matrix_ID")

mean_combined <- bind_rows(mean_list, .id = "Matrix_ID")


# Select relevant columns and convert data from wide to long format
ind_long_format <- ind_combined %>%
                   pivot_longer(cols = -c(Matrix_ID, Sample_ID_1, Population_1, CHRType, K),
                   names_to = "Var2",
                   values_to = "value")

mean_long_format <- mean_combined %>%
                    pivot_longer(cols = -c(Matrix_ID, Population_1, CHRType, K),
                    names_to = "Var2",
                    values_to = "value")


# Renames some columns and select the relevant ones ~ 
fulldf_ind <- ind_long_format %>%
              rename_with(~"Sample_ID_2", Var2) %>%
              rename_with(~"Value", value) %>%
              select(Sample_ID_1, Sample_ID_2, Population_1, CHRType, K, Value) %>%
              ungroup()

fulldf_mean <- mean_long_format %>%
               rename_with(~"Population_2", Var2) %>%
               rename_with(~"Value", value) %>%
               select(Population_1, Population_2, CHRType, K, Value) %>%
               ungroup()


# Gets Population_2 ~
fulldf_ind <- fulldf_ind %>%
              mutate(Population_2 = ifelse(grepl("FR0", Sample_ID_2), "Sales",
                     ifelse(grepl("KAZ", Sample_ID_2), "Chokpak",
                     ifelse(grepl("Lesina", Sample_ID_2), "Lesina",
                     ifelse(grepl("Crotone", Sample_ID_2), "Crotone",
                     ifelse(grepl("Guglionesi", Sample_ID_2), "Guglionesi",
                     ifelse(grepl("PI22NLD0001M", Sample_ID_2), "Y150239",
                     ifelse(grepl("PD22NLD0146F", Sample_ID_2), "Garderen",
                     ifelse(grepl("PD22NLD0147F", Sample_ID_2), "Garderen",
                     ifelse(grepl("PDOM2022NLD0077M", Sample_ID_2), "Meerkerk",
                     ifelse(grepl("PDOM2022NLD0", Sample_ID_2), "Utrecht", "Error"))))))))))) %>%
             select(1:3, Population_2, everything())
  

# Gets Indexes_1 ~
Indexes_1 <- fulldf_ind %>%
  filter(CHRType == "Allosome" & K == "K = 2") %>%
  group_by(Population_1, Sample_ID_1) %>%
  mutate(first_occurrence = row_number() == 1) %>%
  filter(first_occurrence == TRUE) %>%
  group_by(Population_1) %>%
  mutate(Ind_1 = paste0(Population_1, "_", sprintf("%02d", ave(seq_along(Population_1), Population_1, FUN = seq_along)))) %>%
  select(Sample_ID_1, Ind_1) %>%
  ungroup()
Indexes_1 <- Indexes_1 %>%
            select(-Population_1)

# Gets Ind_1 ~
fulldf_ind <- merge(fulldf_ind, Indexes_1, by.x = "Sample_ID_1" , all.x = TRUE)


# Gets Indexes_2 ~
Indexes_2 <- fulldf_ind %>%
  filter(CHRType == "Allosome" & K == "K = 2") %>%
  group_by(Population_2, Sample_ID_2) %>%
  mutate(first_occurrence = row_number() == 1) %>%
  filter(first_occurrence == TRUE) %>%
  group_by(Population_2) %>%
  mutate(Ind_2 = paste0(Population_2, "_", sprintf("%02d", ave(seq_along(Population_2), Population_2, FUN = seq_along)))) %>%
  select(Sample_ID_2, Ind_2) %>%
  ungroup()
Indexes_2 <- Indexes_2 %>%
  select(-Population_2)


# Gets Ind_2 ~
fulldf_ind <- merge(fulldf_ind, Indexes_2, by.x = "Sample_ID_2" , all.x = TRUE)


# Reorders columns ~
fulldf_ind <- fulldf_ind %>%
              select(Sample_ID_1, Sample_ID_2, Population_1, Population_2, Ind_1, Ind_2, CHRType, K, Value)


# Converts NAs into 10s ~
fulldf_ind <- replace(fulldf_ind, is.na(fulldf_ind), 10)


# Filters df for unique combinations and not self-combinations ~
fulldf_ind <- fulldf_ind %>%
              group_by(CHRType, K) %>%
              mutate(combination_sorted = paste(sort(c(Ind_1, Ind_2)), collapse = "_")) %>%
              filter(Ind_1 <= Ind_2) %>%
              select(-combination_sorted) %>%
              #filter(Ind_1 != Ind_2) %>%
              ungroup()

fulldf_mean <- fulldf_mean %>%
               group_by(CHRType, K) %>%
               mutate(combination_sorted = paste(sort(c(Population_1, Population_2)), collapse = "_")) %>%
               filter(Population_1 <= Population_2) %>%
               select(-combination_sorted) %>%
               #filter(Population_1 != Population_2) %>%
               ungroup()

              
# Corrects Population names ~
levels(fulldf_ind$Ind_1 <- sub("Y150239_01", "Y150239", fulldf_ind$Ind_1))
levels(fulldf_ind$Ind_2 <- sub("Y150239_01", "Y150239", fulldf_ind$Ind_2))


# Reorders df based on Populations ~
fulldf_ind <- fulldf_ind %>%
              group_by(CHRType, K) %>%
              arrange(Population_1, Population_2) %>%
              ungroup()


# Combine the data frames
fulldf <- bind_rows(mutate(fulldf_ind, Triangle = "Individual"),
                    mutate(fulldf_mean, Triangle = "Population"))


# Reorders CHRType ~
fulldf$CHRType <- factor(fulldf$CHRType, ordered = T,
                         levels = c("Autosomes",
                                    "Allosome"))


# Reorders K ~
fulldf$K <- factor(fulldf$K, ordered = T,
                   levels = c("K = 7",
                              "K = 6", 
                              "K = 5",
                              "K = 4",
                              "K = 3",
                              "K = 2"))


# Define color palette and breaks
color_palette <- c("#001260", "#EAEDE9", "#601200")
nHalf <- 4
Min <- -.1
Max <- .1
Thresh <- 0

rc1 <- colorRampPalette(colors = color_palette[1:2], space = "Lab")(nHalf)
rc2 <- colorRampPalette(colors = color_palette[2:3], space = "Lab")(nHalf)
rampcols <- c(rc1, rc2)
rampcols[c(nHalf, nHalf+1)] <- rgb(t(col2rgb(color_palette[2])), maxColorValue = 256) 

rb1 <- seq(Min, Thresh, length.out = nHalf + 1)
rb2 <- seq(Thresh, Max, length.out = nHalf + 1)[-1]
rampbreaks <- c(rb1, rb2)


# Creates top plot ~
Ind_Plot <-
  ggplot() + 
  geom_tile(data = subset(fulldf, Triangle == "Individual"), aes(Ind_1, Ind_2, fill = as.numeric(Value)), colour = "#000000")  +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradientn(colors = rampcols, breaks = rampbreaks, limits = c(-.1, .1)) +
  facet_grid(K ~ CHRType, scales = "free", space = "free") +
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
        axis.text.x = element_text(color = "#000000", size = 16, face = "bold", angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(color = "#000000", size = 16, face = "bold"),
        axis.ticks = element_line(color = "#000000", linewidth = .3),
        strip.text = element_text(colour = "#000000", size = 24, face = "bold", family = "Optima"),
        strip.background = element_rect(colour = "#000000", fill = "#d6d6d6", linewidth = .3),
        axis.line = element_line(colour = "#000000", linewidth = .3)) +
  guides(fill = guide_legend(title = "", title.theme = element_text(size = 16, face = "bold"),
                             label.theme = element_text(size = 15), reverse = TRUE))


# Saves plot (Boxplot) ~
ggsave(Ind_Plot, file = "Ind.pdf",
       device = cairo_pdf, limitsize = FALSE, scale = 1, width = 30, height = 50, dpi = 600)


# Creates top plot ~
Mean_Plot <-
  ggplot() + 
  geom_tile(data = subset(fulldf, Triangle == "Population"), aes(Population_1, Population_2, fill = as.numeric(Value)), colour = "#000000")  +
  scale_x_discrete(limits = rev, expand = c(0, 0)) +
  scale_y_discrete(limits = rev, expand = c(0, 0)) +
  scale_fill_gradientn(colors = rampcols, breaks = rampbreaks, limits = c(-.1, .1)) +
  facet_grid(K ~ CHRType, scales = "free", space = "free") +
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
        axis.text.x = element_text(color = "#000000", size = 16, face = "bold", angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(color = "#000000", size = 16, face = "bold"),
        axis.ticks = element_line(color = "#000000", linewidth = .3),
        strip.text = element_text(colour = "#000000", size = 24, face = "bold", family = "Optima"),
        strip.background = element_rect(colour = "#000000", fill = "#d6d6d6", linewidth = .3),
        axis.line = element_line(colour = "#000000", linewidth = .3)) +
  guides(fill = guide_legend(title = "", title.theme = element_text(size = 16, face = "bold"),
                             label.theme = element_text(size = 15), reverse = TRUE))

# Saves plot (Boxplot) ~
ggsave(Mean_Plot, file = "Mean.pdf",
       device = cairo_pdf, limitsize = FALSE, scale = 1, width = 30, height = 50, dpi = 600)
ggsave(CareBears, file = "CareBears.png",
       limitsize = FALSE, scale = 1, width = 30, height = 50, dpi = 600)


Layka <- Ind_Plot + annotation_custom(ggplotGrob(Mean_Plot), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)


# Saves plot (Boxplot) ~
ggsave(Layka, file = "Layka.pdf",
       device = cairo_pdf, limitsize = FALSE, scale = 1, width = 30, height = 50, dpi = 600)


# Increments bottom plot ~
CareBears +
  ggplot(fulldf_mean, aes(Population_1, Population_2, fill = as.numeric(Value))) + 
  geom_tile(colour = "#000000") +
  scale_fill_gradientn(colors = rampcols, breaks = rampbreaks, limits = c(-.1, .1)) +
  geom_tile(data = subset(fulldf_ind, Value == 10), aes(fill = "black")) +
  scale_x_discrete(limits = rev, expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  facet_grid(K ~ CHRType, scales = "free", space = "free") +
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
        strip.text = element_text(colour = "#000000", size = 24, face = "bold", family = "Optima"),
        strip.background = element_rect(colour = "#000000", fill = "#d6d6d6", linewidth = .3),
        axis.line = element_line(colour = "#000000", linewidth = .3)) +
  guides(fill = guide_legend(title = "", title.theme = element_text(size = 16, face = "bold"),
                             label.theme = element_text(size = 15), reverse = TRUE))
              

Corres_2 <- as.matrix(read.table("AllSamples_bcftools.raw.vcf.Filtered.Autosomes.NoKinship.NoTreeSparrow.MAFfiltered.Pruned.K2.corres"))

# Reads the annotation file ~
ids <- read.table("AllSamples_bcftools.raw.vcf.Filtered.Autosomes.NoKinship.NoTreeSparrow.MAFfiltered.Pruned.K2.labels",
                  stringsAsFactors = FALSE, sep = "\t", header = FALSE)

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


plotCorRes(cor_mat = Corres_2, pop = as.vector(ids[, 2]), ord = NULL,
           title = "Evaluation of 1000G admixture proportions with K = 2", max_z = .1, min_z = -.1)


#
##
### The END ~~~~~