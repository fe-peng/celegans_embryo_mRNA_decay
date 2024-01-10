library(tidyverse)
library(ggplot2)
library(broom)
library(data.table)
library(purrr)

# Calculate mRNA half-lives within an individual single-cell RNA-seq biological replicate

## Import the SoupX-corrected total expression counts tsv files for each sample
sample_0M <- read.table("soupx_pilot_sample_0M_totalCounts.txt", sep = "\t", header = TRUE)
sample_20A <- read.table("soupx_pilot_sample_20A_totalCounts.txt", sep = "\t", header = TRUE)
sample_40A <- read.table("soupx_pilot_sample_40A_totalCounts.txt", sep = "\t", header = TRUE)

## Merge dfs together by gene
counts_actd <- merge(sample_0M, sample_20A, by = "gene")
counts_actd <- merge(counts_actd, sample_40A, by = "gene")

## Rename columns in counts_actd so that it's easier to reshape data later on
names(counts_actd)[names(counts_actd) == 'soupx_sample_0M'] <- 'P_0'
names(counts_actd)[names(counts_actd) == 'soupx_sample_20A'] <- 'P_20'
names(counts_actd)[names(counts_actd) == 'soupx_sample_40A'] <- 'P_40'

## Turn first column into row names
counts_actd <- counts_actd %>% remove_rownames %>% tibble::column_to_rownames(var = "gene")

## Only include rows that correspond to ribosomal protein genes (rpl and rps genes)
## Get rid of rows corresponding to mprl, mrps, hrpl, trpl, nrps genes
ribosomal_counts_temp <- counts_actd[rownames(counts_actd)[grep("^mrpl", rownames(counts_actd), invert = TRUE)],]
ribosomal_counts_temp <- ribosomal_counts_temp[rownames(ribosomal_counts_temp)[grep("^mrps", rownames(ribosomal_counts_temp), invert = TRUE)],]
ribosomal_counts_temp <- ribosomal_counts_temp[rownames(ribosomal_counts_temp)[grep("^hrpl", rownames(ribosomal_counts_temp), invert = TRUE)],]
ribosomal_counts_temp <- ribosomal_counts_temp[rownames(ribosomal_counts_temp)[grep("^trpl", rownames(ribosomal_counts_temp), invert = TRUE)],]
ribosomal_counts_temp <- ribosomal_counts_temp[rownames(ribosomal_counts_temp)[grep("^nrps", rownames(ribosomal_counts_temp), invert = TRUE)],]
## Only include rows that correspond to rpl genes
rpl_counts <- ribosomal_counts_temp[rownames(ribosomal_counts_temp)[grep("^rpl", rownames(ribosomal_counts_temp))],]
## Only include rows that correspond to rps genes
rps_counts <- ribosomal_counts_temp[rownames(ribosomal_counts_temp)[grep("^rps", rownames(ribosomal_counts_temp))],]
## Combined rpl and rps counts
ribosomal_counts <- rbind(rpl_counts, rps_counts)

## Filter genes by requiring > 30 UMI in timept 0
gene_counts <- counts_actd
gene_counts_filtered <- gene_counts %>%
  filter(P_0 > 30)

## Create a new row in ribosomal_counts that sums up each column
ribosomal_counts["total" ,] <- colSums(ribosomal_counts)
## Take row with sum and save it as a vector
ribosomal_sum <- as.numeric(as.vector(ribosomal_counts[76,]))

## Divide each gene count by ribosomal sum to normalize the counts between samples
gene_counts_matrix <- data.matrix(gene_counts_filtered, rownames.force = TRUE)
geneCounts_norm <- sweep(gene_counts_matrix, 2, ribosomal_sum, `/`)
## Turn back into dataframe
geneCounts_norm <- as.data.frame(geneCounts_norm)

## Duplicate 0min time point columns
geneCounts_norm$zero = geneCounts_norm$P_0

## Calculate fold change for each time point relative to 0min time point
geneCounts_norm_FC <- geneCounts_norm[ , 1:3]/geneCounts_norm[ , 4]

## Get correction factor for each time point--parameter here is using median ribosomal half-life, as calculated from bulk data
c0 = 1
c20 = 0.5^(20/295.5)
c40 = 0.5^(40/295.5)

## Apply correction factor
geneCounts_norm_FC_corrected <- geneCounts_norm_FC
geneCounts_norm_FC_corrected$P_20 <- geneCounts_norm_FC_corrected$P_20 * c20
geneCounts_norm_FC_corrected$P_40 <- geneCounts_norm_FC_corrected$P_40 * c40

## New row names into it's own column
geneCounts_norm_FC_corrected <- setDT(geneCounts_norm_FC_corrected, keep.rownames = "gene")

## Save dataframe as csv file
write.csv(geneCounts_norm_FC_corrected, "pilot_FC_umi30_medianCorrected.csv", row.names = F)

## Time to reshape the data so it goes from wide to long form
## Then make two new columns for replicate and time point by splitting rep_timept in two
geneCounts_norm_FC_corrected_long <- gather(geneCounts_norm_FC_corrected, key = rep_timept, value = count, -gene) %>%
  separate(rep_timept, into = c("rep", "timept"))

## Make timept column values numeric
geneCounts_norm_FC_corrected_long$timept <- as.numeric(as.character(geneCounts_norm_FC_corrected_long$timept))

## Fit transcript abundance over transcription inhibition time course to exponential decay model
fit_nlsModel <- function(df) nls(count ~ C*exp(-a*timept), data = df, start = list(C = 1, a = 0.01), control = list(maxiter = 100))

fitted <- geneCounts_norm_FC_corrected_long %>% 
  group_nest(gene) %>% # nest by gene
  mutate(fit = map(data, possibly(fit_nlsModel, otherwise = NULL)), tidied = map(fit, tidy), augmented = map(fit, augment))

## Get dataframe with estimates for a and C parameters
param_aC <- fitted %>% 
  unnest(tidied) %>%
  dplyr::select(gene, term, estimate) %>% 
  spread(term, estimate)

## Extract residuals
param_RSS <- fitted %>%
  unnest(augmented) %>%
  dplyr::select(gene, timept, .resid) %>%
  spread(timept, .resid)

## Take square of residuals
param_RSS[,2:4] <- (param_RSS[,2:4])^2

## New column that is sum of residuals
param_RSS$RSS <- rowSums(param_RSS[ , 2:4])

## Subset so that only gene name and RSS columns are kept for eventual merge of dataframes
param_RSS2 <- param_RSS[, c("gene", "RSS")]

## Extract counts, ie the actual values
param_TSS <- fitted %>%
  unnest(augmented) %>%
  dplyr::select(gene, timept, count) %>%
  spread(timept, count)

## New column that is the mean of all count values per gene
param_TSS$mean_y <- rowMeans(param_TSS[,-1])

## Substract mean_y from each count and square it
param_TSS[,2] <- ((param_TSS[,2])-(param_TSS[,5]))^2
param_TSS[,3] <- ((param_TSS[,3])-(param_TSS[,5]))^2
param_TSS[,4] <- ((param_TSS[,4])-(param_TSS[,5]))^2

## Get rid of mean_y column...
param_TSS$mean_y <- NULL

## New column TSS that is sum of squared values
param_TSS$TSS <- rowSums(param_TSS[ ,-1])

## Subset so that only gene name and TSS columns are kept for eventual merge of dataframes
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
  mutate(ribosomal_corrected_half_life = (log(2))/a)

write.table(param_all, file = "pilotRep_umi30_correctRibosomalNormalizedHalflife.txt", sep = "\t", row.names = FALSE, quote = FALSE)
