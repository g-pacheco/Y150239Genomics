# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load required packages
pacman::p_load(tidyverse)

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

# Define a function to compute individual correlations
compute_individual_correlations <- function(corres_ind_list, pop_list) {
  result <- lapply(1:length(corres_ind_list), function(i) {
    corres_ind <- corres_ind_list[[i]]
    pop <- pop_list[[i]]
    
    # Initialize a matrix to store individual correlations
    indiv_cors <- matrix(NA, nrow = length(pop), ncol = length(pop))
    colnames(indiv_cors) <- pop
    rownames(indiv_cors) <- pop
    
    # Loop through each pair of populations
    for (p1 in unique(pop)) {
      for (p2 in unique(pop)) {
        # Extract correlations corresponding to the current pair of populations
        cor_values <- corres_ind[pop == p1, pop == p2]
        
        # Compute individual correlation, ignoring NAs
        indiv_cor <- mean(cor_values, na.rm = TRUE)
        
        # Store the individual correlation in the matrix
        indiv_cors[p1, p2] <- indiv_cor
      }
    }
    
    # Return the matrix of individual correlations
    return(indiv_cors)
  })
  
  return(result)
}

# Apply the function to compute individual correlations
indiv_cor_list <- compute_individual_correlations(corres_ind_list = corres, pop_list = lapply(annot, function(ann) ann$Population_1))


exclude_cols <- c("Sample_ID_1", "Population_1", "CHRType", "K")

compute_mean_correlation <- function(corres_ind_list, exclude_cols) {
  result <- lapply(corres_ind_list, function(corres_ind) {
    pop <- corres_ind$Population_1
    
    # Exclude specified columns
    if (!is.null(exclude_cols)) {
      corres_sub <- corres_ind[, !(names(corres_ind) %in% exclude_cols)]
    } else {
      corres_sub <- corres_ind
    }
    
    # Filter out non-numeric columns
    numeric_cols <- sapply(corres_sub, is.numeric)
    corres_sub <- corres_sub[, numeric_cols]
    
    mean_cors <- matrix(NA, ncol = length(unique(pop)), nrow = length(unique(pop)))
    colnames(mean_cors) <- unique(pop)
    rownames(mean_cors) <- unique(pop)
    
    for (i1 in 1:length(unique(pop))) {
      for (i2 in 1:length(unique(pop))) {
        p1 <- unique(pop)[i1]
        p2 <- unique(pop)[i2]
        
        # Filter out rows corresponding to populations p1 and p2
        corres_sub_filtered <- corres_sub[pop == p1, pop == p2]
        
        print("Correspondence matrix for populations:")
        print(corres_sub_filtered)
        
        # Compute mean if there are enough non-NA values
        if (sum(!is.na(corres_sub_filtered)) > 1) {
          mean_cors[i1, i2] <- mean(corres_sub_filtered, na.rm = TRUE)
        } else {
          mean_cors[i1, i2] <- 0
        }
      }
    }
    return(mean_cors)
  })
  return(result)
}





mean_cors_list <- compute_mean_correlation(main_list,
                                           exclude_cols = c("Sample_ID_1", "Population_1", "CHRType", "K"))

# Function to compute individual correlations
compute_on_matrices <- function(corres_ind_list, pop_list, ord_list) {
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
    indiv_matrix <- cbind(annotation_df, corres_ind)
    return(indiv_matrix)})
  return(result)}

# Apply functions
main_list <- compute_on_matrices(corres_ind_list = corres,
                                 pop_list = lapply(annot, function(ann) ann$Population_1),
                                 ord_list = lapply(corres, orderInds))

# Convert all matrices to numeric type
main_list <- lapply(main_list, as.matrix)
main_list <- lapply(main_list, as.data.frame)

# Combine all matrices into a single dataset
combined_data <- bind_rows(main_list, .id = "Matrix_ID")

# Select relevant columns and convert data from wide to long format
long_format <- combined_data %>%
  pivot_longer(cols = -c(Matrix_ID, Sample_ID_1, Population_1, CHRType, K),
               names_to = "Var2",
               values_to = "value") %>%
  rename_with(~ "Sample_ID_2", Var2) %>%
  rename_with(~ "Value", value) %>%
  select(Sample_ID_1, Sample_ID_2, Population_1, CHRType, K, Value) %>%
  ungroup()

# Get Population_2
fulldf_ind <- long_format %>%
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

# Get Indexes_1
Indexes_1 <- fulldf_ind %>%
  filter(CHRType == "Allosome" & K == "K = 2") %>%
  group_by(Population_1, Sample_ID_1) %>%
  mutate(first_occurrence = row_number() == 1) %>%
  filter(first_occurrence == TRUE) %>%
  group_by(Population_1) %>%
  mutate(Ind_1 = paste0(Population_1, "_", sprintf("%02d", ave(seq_along(Population_1), Population_1, FUN = seq_along)))) %>%
  select(Sample_ID_1, Ind_1) %>%
  ungroup() %>%
  select(-Population_1)

# Get Ind_1
fulldf_ind <- merge(fulldf_ind, Indexes_1, by.x = "Sample_ID_1", all.x = TRUE)

# Get Indexes_2
Indexes_2 <- fulldf_ind %>%
  filter(CHRType == "Allosome" & K == "K = 2") %>%
  group_by(Population_2, Sample_ID_2) %>%
  mutate(first_occurrence = row_number() == 1) %>%
  filter(first_occurrence == TRUE) %>%
  group_by(Population_2
           