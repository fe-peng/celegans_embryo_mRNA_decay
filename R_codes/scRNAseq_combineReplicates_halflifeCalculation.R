library(tidyverse)
library(ggplot2)
library(broom)
library(data.table)
library(purrr)

# Calculate mRNA half-lives across 3 biological replicates

## Import normalized UMI for each replicate
pilot_FC <- read.csv(file = 'pilot_final/pilot_FC_umi30_medianCorrected.csv')
second_FC <- read.csv(file = '12.10_rep_final/second_FC_umi30_medianCorrected.csv')
## Rename columns 
names(second_FC)[names(second_FC) == 'P_0'] <- 'S_0'
names(second_FC)[names(second_FC) == 'P_20'] <- 'S_20'
names(second_FC)[names(second_FC) == 'P_40'] <- 'S_40'
third_FC <- read.csv(file = '12.14_rep_final/third_FC_umi30_medianCorrected.csv')
## Rename columns 
names(third_FC)[names(third_FC) == 'P_0'] <- 'T_0'
names(third_FC)[names(third_FC) == 'P_20'] <- 'T_20'
names(third_FC)[names(third_FC) == 'P_40'] <- 'T_40'

## Merge reps to only keep genes shared in common
all_FC <- merge(pilot_FC, second_FC, by = "gene")
all_FC <- merge(all_FC, third_FC, by = "gene")

## Filter individual replicate FC dataframes to only keep genes kept between all 3
pilot_FCshared <- pilot_FC %>%
  filter(gene %in% all_FC$gene)

second_FCshared <- second_FC %>%
  filter(gene %in% all_FC$gene)

third_FCshared <- third_FC %>%
  filter(gene %in% all_FC$gene)

## Time to reshape the data so it goes from wide to long form
## Then make two new columns for replicate and time point by splitting rep_timept in two
pilot_FCshared_long <- gather(pilot_FCshared, key = rep_timept, value = count, -gene) %>%
  separate(rep_timept, into = c("rep", "timept"))

second_FCshared_long <- gather(second_FCshared, key = rep_timept, value = count, -gene) %>%
  separate(rep_timept, into = c("rep", "timept"))

third_FCshared_long <- gather(third_FCshared, key = rep_timept, value = count, -gene) %>%
  separate(rep_timept, into = c("rep", "timept"))

## Combine dataframes vertically
all_FC_long <- rbind(pilot_FCshared_long, second_FCshared_long)
all_FC_long <- rbind(all_FC_long, third_FCshared_long)

## Make timept column values numeric
all_FC_long$timept <- as.numeric(as.character(all_FC_long$timept))

## Fit transcript abundance over transcription inhibition time course to exponential decay model
fit_nlsModel <- function(df) nls(count ~ C*exp(-a*timept), data = df, start = list(C = 1, a = 0.01), control = list(maxiter = 100))

fitted <- all_FC_long %>% 
  group_nest(gene) %>% # nest by gene
  mutate(fit = map(data, possibly(fit_nlsModel, otherwise = NULL)), tidied = map(fit, tidy), augmented = map(fit, augment))

## Get dataframe with estimates for a and C parameters
param_aC <- fitted %>% 
  unnest(tidied) %>%
  dplyr::select(gene, term, estimate) %>% 
  spread(term, estimate)

## Extract residuals
param_RSS <- fitted %>%
  unnest(augmented)

## Create extra column where each row for each gene gets a different number
param_RSS <- param_RSS %>%
  group_by(gene) %>%
  mutate(grouped_id = row_number())

param_RSS <- param_RSS %>%
  dplyr::select(gene, grouped_id, .resid) %>%
  spread(grouped_id, .resid)

## Take square of residuals
param_RSS[,2:10] <- (param_RSS[,2:10])^2

## New column that is sum of residuals
param_RSS$RSS <- rowSums(param_RSS[ , 2:10])

## Subset so that only gene name and RSS columns are kept for eventual merge of dataframes
param_RSS2 <- param_RSS[, c("gene", "RSS")]

## Extract counts, ie the actual values
## Unnest augmented column
param_TSS <- fitted %>%
  unnest(augmented) 

## Create extra column where each row for each gene gets a different number
param_TSS <- param_TSS %>%
  group_by(gene) %>%
  mutate(grouped_id = row_number())

param_TSS <- param_TSS %>%
  dplyr::select(gene, grouped_id, count) %>%
  spread(grouped_id, count)

## New column that is the mean of all count values per gene
param_TSS$mean_y <- rowMeans(param_TSS[,-1])

## Substract mean_y from each count and square it...
param_TSS[,2] <- ((param_TSS[,2])-(param_TSS[,11]))^2
param_TSS[,3] <- ((param_TSS[,3])-(param_TSS[,11]))^2
param_TSS[,4] <- ((param_TSS[,4])-(param_TSS[,11]))^2
param_TSS[,5] <- ((param_TSS[,5])-(param_TSS[,11]))^2
param_TSS[,6] <- ((param_TSS[,6])-(param_TSS[,11]))^2
param_TSS[,7] <- ((param_TSS[,7])-(param_TSS[,11]))^2
param_TSS[,8] <- ((param_TSS[,8])-(param_TSS[,11]))^2
param_TSS[,9] <- ((param_TSS[,9])-(param_TSS[,11]))^2
param_TSS[,10] <- ((param_TSS[,10])-(param_TSS[,11]))^2

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

## Exclude the few genes where 'a' is negative
param_all <- subset(param_all, a > 0)

## Calculate half-lives and put into new column
param_all <- param_all %>% 
  mutate(half_life = (log(2))/a)

## Calculate coefficient of variation for half-lives
## Import half-lives calculated for individual bio replicates
pilot_hl <- read.table("pilotRep_umi30_correctRibosomalNormalizedHalflife.txt", sep = "\t", header = TRUE)
second_hl <- read.table("secondRep_umi30_correctRibosomalNormalizedHalflife.txt", sep = "\t", header = TRUE)
third_hl <- read.table("thirdRep_umi30_correctRibosomalNormalizedHalflife.txt", sep = "\t", header = TRUE)

## Only keep 'gene' and 'half_life' columns
pilot_hl <- pilot_hl[ , c("gene", "ribosomal_corrected_half_life")] 
second_hl <- second_hl[ , c("gene", "ribosomal_corrected_half_life")] 
third_hl <- third_hl[ , c("gene", "ribosomal_corrected_half_life")] 

## Rename columns so that half-lives reflect which rep they came from 
names(pilot_hl)[names(pilot_hl) == 'ribosomal_corrected_half_life'] <- 'ribosomal_corrected_half_life_pilot'
names(second_hl)[names(second_hl) == 'ribosomal_corrected_half_life'] <- 'ribosomal_corrected_half_life_second'
names(third_hl)[names(third_hl) == 'ribosomal_corrected_half_life'] <- 'ribosomal_corrected_half_life_third'

## Merge all dataframes by gene to only keep genes shared across all reps
all_CV <- merge(pilot_hl, second_hl, by = "gene")
all_CV <- merge(all_CV, third_hl, by = "gene")

## Calculate standard deviation and mean for each gene across replicates ABD
all_CV <- transform(all_CV, stdev = apply(all_CV[, 2:4], 1, sd))
all_CV <- transform(all_CV, mean = apply(all_CV[, 2:4], 1, mean))

## Divide stdev by mean half-life to get a measure of reproducibility of half-lives between replicates (coefficient of variation)
all_CV <- transform(all_CV, CV = stdev / mean)

library(nlstools)

## Turn param_all into a tibble so it can be merged with fitted, a tibble
param_all_tbl <- as_tibble(param_all)

## merge the two tables
CI_test <- merge(param_all_tbl, fitted, by = "gene")

## Only keep 'gene' and 'fit' columns
CI_test <- CI_test[ , c("gene", "fit")]

## Use nlstools::confint2 to calculate confidence intervals
CI_test <- CI_test %>%
  mutate(confidence_interval = map(fit, nlstools::confint2, level = 0.95))

## Create new column 'transposed' that is the CI matrices transposed
CI_test <- CI_test %>%
  mutate(transposed = map(confidence_interval, t))

## Turn each matrix into a vector
CI_test <- CI_test %>%
  mutate(transposed_rev = map(transposed, rev))

## Create new dataframe that is just the 'gene' and 'transposed_rev' columns in CI_test
confidence_intervals95 <- CI_test[ , c("gene", "transposed_rev")]

## Basically, take the four different values in 'transposed_rev' and turn each one into its own column 
confidence_intervals95 <- confidence_intervals95 %>%
  unnest(transposed_rev) %>%
  group_by(gene) %>%
  mutate(key = row_number()) %>% 
  spread(key, transposed_rev)

## Rename the columns so they indicate which variable and which CI limit they correspond to
names(confidence_intervals95)[names(confidence_intervals95) == '1'] <- 'a_97.5'
names(confidence_intervals95)[names(confidence_intervals95) == '2'] <- 'a_2.5'
names(confidence_intervals95)[names(confidence_intervals95) == '3'] <- 'C_97.5'
names(confidence_intervals95)[names(confidence_intervals95) == '4'] <- 'C_2.5'

## Create new columns that translate the 95% CI for the decay rate into 95% CI for half-life values
confidence_intervals95 <- confidence_intervals95 %>%
  mutate(half_life_97.5 = (log(2))/a_97.5)

confidence_intervals95 <- confidence_intervals95 %>%
  mutate(half_life_2.5 = (log(2))/a_2.5)

## Merge dataframes together so that get one big dataframe with variable estimates, R^2, half-live values, and more
param_all_CI <- merge(param_all, confidence_intervals95, by = "gene")

## Create new column FC that is half_life_2.5/half_life
param_all_CI$FC_2.5 <- param_all_CI$half_life_2.5 / param_all_CI$half_life

## Create new column FC that is half_life/half_life_97.5
param_all_CI$FC_97.5 <- param_all_CI$half_life / param_all_CI$half_life_97.5 

## Subset
param_all_CI <- param_all_CI[ , c("gene", "half_life", "half_life_97.5", "half_life_2.5", "rsq", "FC_2.5", "FC_97.5")]

## Merge to get CI and CV data into same df
CI_CV <- merge(param_all_CI, all_CV, by = "gene")

## Implement moderate filtering strategy for half-lives
## Confidence interval filtering
CI_FC3 <- subset(CI_CV, FC_2.5 >= 0 & FC_2.5 <= 3)
CI_FC3 <- subset(CI_FC3, half_life >= 0)
## Coefficient of variation filtering
CV50 <- subset(CI_CV, CV <= 0.5)
CV50 <- subset(CV50, CV >= 0)
## CI or CV moderate filtering
CIFC3_CV50 <- rbind(CI_FC3, CV50)
CIFC3_CV50 <- unique(CIFC3_CV50)

## Looser filtering strategy for genes with long half-lives
CV75 <- subset(CI_CV, CV <= 0.75)
CV75 <- subset(CV75, CV >= 0)
CI_CV_100min <- subset(CI_CV, half_life > 100)
CV75_hl100 <- subset(CV75, WB_ID %in% CI_CV_100min$WB_ID)
CI_FC4 <- subset(CI_CV, FC_2.5 >= 0 & FC_2.5 <= 4)
CI_FC4 <- subset(CI_FC4, half_life >= 0)
CIFC4_hl100 <- subset(CI_CV_100min, WB_ID %in% CI_FC4$WB_ID)

## Get a list of genes with CV50 or CIFC3, and also half-lives > 100 min with CV75 or CIFC4 threshold
CIFC3_CV50_CIFC4andCV75hl100 <- rbind(CIFC3_CV50, CV75_hl100, CIFC4_hl100)
CIFC3_CV50_CIFC4andCV75hl100 <- unique(CIFC3_CV50_CIFC4andCV75hl100)

write.table(CIFC3_CV50_CIFC4andCV75hl100, file = "pseudobulk_halflives_CIFC3_CV50_CIFC4andCV75hl100_umi30.txt", sep = "\t", row.names = TRUE, quote = FALSE)
