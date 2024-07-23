### The BEGINNING ~~~~~
##
# ~ Plots Runs visFuns.R | Provided by evalAdmix.R.


# Cleans the environment ~ 
rm(list=ls())


# Sets working directory ~
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Loads required packages ~
pacman::p_load(tidyverse, ggnewscale)
source("visFuns.R")

# Reads data & annotation ~
Corres_2 <- as.matrix(read.table("AllSamples_bcftools.raw.vcf.Filtered.Allosome.NoKinship.NoTreeSparrow.MAFfiltered.Pruned.K2.corres"))

ids <- read.table("AllSamples_bcftools.raw.vcf.Filtered.Allosome.NoKinship.NoTreeSparrow.MAFfiltered.Pruned.K2.labels",
                  stringsAsFactors = FALSE, sep = "\t", header = FALSE)

q <- read.table("AllSamples_bcftools.raw.vcf.Filtered.Allosome.NoKinship.NoTreeSparrow.MAFfiltered.Pruned.K2.qopt",
                stringsAsFactors = TRUE)


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


#ord_without_q <- orderInds(pop = as.vector(ids[, 2]))
ord_with_q <- orderInds(pop = as.vector(ids[, 2]), q = q)

quartz(width = 10, height = 10)
plotCorRes(cor_mat = Corres_2, pop = as.vector(ids[, 2]),
           title = "Evaluation of 1000G admixture proportions with K = 2", max_z = .3, min_z = -.3)


Corres_2_Pop <- as.matrix(read.table("mean_cors.csv", sep = ",", header = TRUE, row.names = 1))
Corres_2_Final <- as.matrix(read.table("cor_mat.csv", sep = ",", header = TRUE, row.names = 1))
