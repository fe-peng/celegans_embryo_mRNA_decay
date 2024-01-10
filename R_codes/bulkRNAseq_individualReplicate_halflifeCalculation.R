library(tidyverse)
library(dplyr)
library(broom)
library(data.table)
library(purrr)

# Calculate mRNA half-lives within an individual bulk RNA-seq biological replicate

## Import counts data for replicate A
repA_counts <- read.csv(file = 'repA_counts.csv')

## Only include rows that are genes, by excluding rows corresponding to ERCC spike-in RNAs
gene_counts <- repA_counts[rownames(repA_counts)[grep("^ERCC", rownames(repA_counts), invert = TRUE)],]

## Only include rows that correspond to ERCC spike-in RNAs
spike_counts <- repA_counts[rownames(repA_counts)[grep("^ERCC", rownames(repA_counts))],]

## Require genes to have count > 30 in time pt 0
gene_counts_filtered <- gene_counts %>%
  filter(A_0 > 30)

## Create a new row in spike_counts that sums up each column
spike_counts["total" ,] <- colSums(spike_counts)

## Take row with sum and save it as a vector
spike_sum <- as.numeric(as.vector(spike_counts[93,]))

## Turn into matrix
gene_counts_matrix <- data.matrix(gene_counts_filtered, rownames.force = TRUE)

## Divide each gene count by spike-in sum to normalize the counts between samples
geneCounts_norm <- sweep(gene_counts_matrix, 2, spike_sum, `/`)

## Turn back into dataframe
geneCounts_norm <- as.data.frame(geneCounts_norm)

## Duplicate time pt 0 column
geneCounts_norm$zero = geneCounts_norm$A_0

## Calculate fold change for each time point relative to time pt 0
geneCounts_norm_FC <- geneCounts_norm[ , 1:5]/geneCounts_norm[ , 6]

## Row names into its own column...
geneCounts_norm_FC <- setDT(geneCounts_norm_FC, keep.rownames = "gene")

## Save dataframe as csv file
write.csv(geneCounts_norm_FC, "repA_counts30_FC.csv", row.names = F)

## Time to reshape the data so it goes from wide to long form
## Then make two new columns for replicate and time point by splitting rep_timept in two
geneCounts_norm_FC_long <- gather(geneCounts_norm_FC, key = rep_timept, value = count, -gene) %>%
  separate(rep_timept, into = c("rep", "timept"))

## Make timept column values numeric
geneCounts_norm_FC_long$timept <- as.numeric(as.character(geneCounts_norm_FC_long$timept))

## Fit transcript abundance over transcription inhibition time course to exponential decay model
fit_nlsModel <- function(df) nls(count ~ C*exp(-a*timept), data = df, start = list(C = 1, a = 0.01), control = list(maxiter = 100))

fitted <- geneCounts_norm_FC_long %>% 
  group_nest(gene) %>% # nest by gene
  mutate(fit = map(data, possibly(fit_nlsModel, otherwise = NULL)), tidied = map(fit, tidy), augmented = map(fit, augment)) 

## Get dataframe with estimates for a and C parameters
param_aC <- fitted %>% 
  unnest(tidied) %>%
  select(gene, term, estimate) %>% 
  spread(term, estimate)

## Extract residuals
param_RSS <- fitted %>%
  unnest(augmented) %>%
  select(gene, timept, .resid) %>%
  spread(timept, .resid)

## Take square of residuals
param_RSS[,2:6] <- (param_RSS[,2:6])^2

## New column that is sum of residuals
param_RSS$RSS <- rowSums(param_RSS[ , 2:6])

## Subset so that only gene name and RSS columns are kept for eventual merge of dataframes
param_RSS2 <- param_RSS[, c("gene", "RSS")]

## Extract counts, ie the actual values
param_TSS <- fitted %>%
  unnest(augmented) %>%
  select(gene, timept, count) %>%
  spread(timept, count)

## New column that is the mean of all count values per gene
param_TSS$mean_y <- rowMeans(param_TSS[,-1])

## Substract mean_y from each count and square it
param_TSS[,2] <- ((param_TSS[,2])-(param_TSS[,7]))^2
param_TSS[,3] <- ((param_TSS[,3])-(param_TSS[,7]))^2
param_TSS[,4] <- ((param_TSS[,4])-(param_TSS[,7]))^2
param_TSS[,5] <- ((param_TSS[,5])-(param_TSS[,7]))^2
param_TSS[,6] <- ((param_TSS[,6])-(param_TSS[,7]))^2

## Get rid of mean_y column...
param_TSS$mean_y <- NULL

## New column TSS that is sum of squared values
param_TSS$TSS <- rowSums(param_TSS[ ,-1])

## subset so that only gene name and TSS columns are kept for eventual merge of dataframes
param_TSS2 <- param_TSS[ , c("gene", "TSS")]

## Merge tables by gene column
param_rsq <- merge(param_RSS2, param_TSS2, by = "gene")

## New column that is 1 - RSS/TSS to get R2 fit value for each gene to exponential decay
param_rsq$rsq <- 1 - (param_rsq[,2]/param_rsq[,3])

## Merge tables by gene column
param_all <- merge(param_aC, param_rsq, by = "gene")

## Exclude the few genes where 'a' is negative...
param_all <- subset(param_all, a > 0)

## Calculate half-lives and put into new column
param_all <- param_all %>% 
  mutate(half_life = (log(2))/a)

## Save dataframe as csv file
write.csv(param_all, "repA_nls_halflives_count30.csv", row.names = F)
