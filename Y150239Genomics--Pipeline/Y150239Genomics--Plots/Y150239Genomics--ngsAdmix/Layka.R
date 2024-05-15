### The BEGINNING ~~~~~
##
# ~ Plots NLSparrow--evalAdmix | First written by Jose Samaniego with later modifications by George Pacheco.


# Cleans the environment ~ 
rm(list=ls())


# Sets working directory ~
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Loads required packages ~
pacman::p_load(tidyverse, reshape2, pheatmap, scales, optparse, plyr, RColorBrewer, extrafont)
source("visFuns.R")


# Imports data while incorporating annotation ~
corres <- list()
annot <- list()
corres_files <- dir(pattern = ".corres")
annot_files <- dir(pattern = ".labels")
for (k in seq_along(annot_files)) {
  annot[[k]] <- read.table(annot_files[k], sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  colnames(annot[[k]]) <- c("Sample_ID")
  annot[[k]]$Population <- ifelse(grepl("FR0", annot[[k]]$Sample_ID), "Sales",
                           ifelse(grepl("KAZ", annot[[k]]$Sample_ID), "Chokpak",
                           ifelse(grepl("Lesina", annot[[k]]$Sample_ID), "Lesina",
                           ifelse(grepl("Crotone", annot[[k]]$Sample_ID), "Crotone",
                           ifelse(grepl("Guglionesi", annot[[k]]$Sample_ID), "Guglionesi",
                           ifelse(grepl("PI22NLD0001M", annot[[k]]$Sample_ID), "Y150239",
                           ifelse(grepl("PD22NLD0146F", annot[[k]]$Sample_ID), "Garderen",
                           ifelse(grepl("PD22NLD0147F", annot[[k]]$Sample_ID), "Garderen",
                           ifelse(grepl("PDOM2022NLD0077M", annot[[k]]$Sample_ID), "Meerkerk",
                           ifelse(grepl("PDOM2022NLD0", annot[[k]]$Sample_ID), "Utrecht", "Error"))))))))))
  corres_df <- as.data.frame(read.table(corres_files[k]))
  rownames(corres_df) <- annot[[k]]$Sample_ID
  colnames(corres_df) <- annot[[k]]$Sample_ID
  corres_df$CHRType <- str_extract(corresL[k], "(Allosome|Autosomes)")
  corres[[k]] <- corres_df}


Oi <- corres[1]


# Performs mean computation ~
compute_mean_correlation <- function(corres_ind, pop) {
  N <- length(pop)
  mean_cors <- matrix(NA, ncol = length(unique(pop)), nrow = length(unique(pop)))
  colnames(mean_cors) <- unique(pop)
  rownames(mean_cors) <- unique(pop)
  for (i1 in 1:(length(unique(pop)))) {
    for (i2 in 1:(length(unique(pop)))) {
      p1 <- unique(pop)[i1]
      p2 <- unique(pop)[i2]
      mean_cors[i1, i2] <- mean(corres_ind[which(pop == p1), which(pop == p2)], na.rm = TRUE)}}
  return(mean_cors)}


# Perform individual computation ~
compute_on_matrices <- function(corres_ind_list, pop_list, ord_list) {
  result <- lapply(1:length(corres_ind_list), function(i) {
    corres_ind <- corres_ind_list[[i]]
    pop <- pop_list[[i]]
    ord <- ord_list[[i]]
    mean_cor <- compute_mean_correlation(corres_ind, pop)
    indiv_cor <- list(
      individual = corres_ind,
      mean_correlation = mean_cor)
    return(indiv_cor)})
  return(result)}


# Applies functions ~
main_list <- compute_on_matrices(corres_ind_list = corres,
                                 pop_list = lapply(annot, function(ann) ann$Population),
                                 ord_list = lapply(corres, orderInds))








# Define the pattern for the condition
pattern <- "individual"

# Define values for the new column based on the condition
value_X <- "Individual"
value_Y <- "Population"

# Loop through each sublist in the main list
for (sublist in main_list) {
  # Loop through each matrix in the sublist
  for (matrix_name in names(sublist)) {
    matrix <- sublist[[matrix_name]]
    # Check if the matrix name matches the pattern
    if (grepl(pattern, matrix_name)) {
      # Add a new column with value X
      matrix$New_Column <- value_X
    } else {
      # Add a new column with value Y
      matrix$New_Column <- value_Y
    }
    # Assign the modified matrix back to the sublist
    sublist[[matrix_name]] <- matrix
  }
}


first_sublist <- main_list[[1]]
first_matrix <- first_sublist[[1]]



# Check if the matrix name matches the pattern
if (grepl(pattern, names(matrix))) {
  # Add new columns with value "X" to the matrix
  matrix[["Computation"]] <- "Individual"
} else {
  # Add new columns with value "Y" to the matrix
  matrix[["Computation"]] <- "Population"}


