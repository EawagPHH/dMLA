#dMLA Bangladesh strains-------------------
##Load packages--------
library (ggplot2)
library(dplyr)
library(gridExtra)
library(ggpubr)
library(tidyr)
library(seqinr)

##Formatting df------
#Import results
setwd("~/Desktop/dmla-amr-vfs/data")
reads <- read.csv("ww_output.csv", header = TRUE)
reads

reads <- reads %>% 
  mutate(Sample_real = as.character(X1)) %>% 
  mutate(Sample_real = case_when(
    Sample_real == "CCACTAG" ~ "48_1",
    Sample_real == "GTCTTCT" ~ "49_1",Sample_real == "ACAAAGC" ~ "50_1",Sample_real == "AAAAGGC" ~ "51_1",Sample_real == "ACTGTGT" ~ "52_1", Sample_real == "CTGGTAC" ~ "53_1",
    Sample_real == "GTTGCTA" ~ "54_1", Sample_real == "CCATTCA" ~ "zNegative_1",Sample_real == "GATATCG" ~ "55_1",Sample_real == "ACAGGAT" ~ "56_1",Sample_real == "CGCATAC" ~ "57_1",
    Sample_real == "GGACCTA" ~ "58_1", Sample_real == "GAACTGA" ~ "59_1", Sample_real == "AGTCGTG" ~ "60_1",Sample_real == "TTCCAGG" ~ "61_1",Sample_real == "TAGAGCG" ~ "zNegative_2",
    Sample_real == "AAACCCT" ~ "62_1",Sample_real == "ATCCAGT" ~ "63_1",Sample_real == "TCTTCGT" ~ "64_1",Sample_real == "GTGTCCT" ~ "65_1",Sample_real == "CACAGAT" ~ "66_1",
    Sample_real == "ATAGAGC" ~ "68_1", Sample_real == "GGTCATG" ~ "69_1", Sample_real == "CCATCTC" ~ "zNegative_3",Sample_real == "TATGCAG" ~ "70_1",Sample_real == "AACCGCT" ~ "71_1",
    Sample_real == "TGACTCT" ~ "72_1",Sample_real == "AAGCACA" ~ "73_1",Sample_real == "TGTGAAC" ~ "75_1",Sample_real == "AACAGTG" ~ "76_1",Sample_real == "GCCTATA" ~ "77_1",
    Sample_real == "GTCCTCT" ~ "zNegative_4",Sample_real == "TAGAACG" ~ "48_2",Sample_real == "TGGGTGA" ~ "49_2", Sample_real == "ACACGTG" ~ "50_2",Sample_real == "TGCCCAA" ~ "51_2",
    Sample_real == "TCGTACA" ~ "52_2", Sample_real == "AACCAAG" ~ "53_2", Sample_real == "GTGCTAT" ~ "54_2",Sample_real == "TCCTCAT" ~ "zNegative_5",Sample_real == "GTCGTAT" ~ "55_2",
    Sample_real == "TCTGGAA" ~ "56_2",Sample_real == "CAGTAGG" ~ "57_2", Sample_real == "ACGTCAT" ~ "58_2",Sample_real == "ACGGAAC" ~ "59_2",Sample_real == "CGGATGA" ~ "60_2",
    Sample_real == "CCGCATA" ~ "61_2", Sample_real == "GTGAAGA" ~ "zNegative_6", Sample_real == "AACGTGT" ~ "62_2",Sample_real == "AGAAGAC" ~ "63_2",Sample_real == "TGTAGGG" ~ "64_2",
    Sample_real == "AATGCCT" ~ "65_2", Sample_real == "GCTTTCT" ~ "66_2", Sample_real == "CATGAAG" ~ "68_2",Sample_real == "CGATCAC" ~ "69_2",Sample_real == "TCGCGTT" ~ "zNegative_7",
    Sample_real == "GATTGGC" ~ "aggR11_1", Sample_real == "CATCGGT" ~ "71_2", Sample_real == "CTTTCCA" ~ "72_2",Sample_real == "ACCAGAT" ~ "73_2",Sample_real == "GTGATAC" ~ "75_2",
    Sample_real == "TGCATCC" ~ "76_2", Sample_real == "TGCAAAC" ~ "77_2", Sample_real == "GCAACGA" ~ "zNegative_8",Sample_real == "GCCATAC" ~ "70_2",Sample_real == "TACCTTC" ~ "aggR11_2",
    Sample_real == "GAATCGA" ~ "aaiC47_1", Sample_real == "CTACGTT" ~ "aaiC47_2", Sample_real == "GTTTCGG" ~ "zPCRNeg_1",Sample_real == "GCAAATG" ~ "zPCRNeg_2")) 

#Distribution of reads------------------
##Filter of negative samples containing positive
# Define the list of negative samples
negative_samples <- c("zNegative_1", "zNegative_2", "zNegative_3", "zNegative_4",
                      "zNegative_5", "zNegative_6", "zNegative_7", "zNegative_8",
                      "zPCRNeg_1", "zPCRNeg_2")

# Filter the dataframe
filtered_reads <- reads %>%
  filter(Sample_real %in% negative_samples |
           (Sample_real %in% c("aggR11_1", "aggR11_2") & X2 != "aggR11") |
           (Sample_real %in% c("aaiC47_1", "aaiC47_2") & X2 != "aaiC47"))

#Plot distribution of false positive for each target
ggplot(filtered_reads, aes(x=n)) + 
  geom_histogram() +
  theme_minimal() +
  labs(title="Distribution of of false positive read counts for each probe-pair - wastewater isolates",
       x="n",
       y="Frequency") +
  theme(legend.title = element_blank())+
 facet_wrap(X2 ~ ., ncol = 10, scales="free")

#Calculate the 99.9% quantile of the gamma distribution of each the false positive for each target as threshold
# Assuming 'filtered_reads' is your dataframe with columns 'n' and 'X2'
results_threshold <- filtered_reads %>%
  group_by(X2) %>%
  summarise(
    mean_data = mean(n),
    var_data = var(n)
  ) %>%
  mutate(
    shape_estimate = mean_data^2 / var_data,
    rate_estimate = mean_data / var_data,
    quantile_99_9 = qgamma(0.999, shape = shape_estimate, rate = rate_estimate)
  )

# Step 1: Calculate a fallback value for each X2 in filtered_reads
fallback_values <- filtered_reads %>%
  group_by(X2) %>%
  summarise(fallback_n = max(n), .groups = 'drop') # Using max, adjust as needed

# Step 2: Join fallback_values with results_threshold
results_threshold <- results_threshold %>%
  left_join(fallback_values, by = "X2")

# Step 3: Replace NA in quantile_99_9 with fallback_n
results_threshold <- results_threshold %>%
  mutate(
    quantile_99_9 = ifelse(is.na(quantile_99_9), fallback_n, quantile_99_9)
  ) %>% 
  dplyr::select(-fallback_n) # Removing fallback_n to clean up the dataframe

#Join results_threshold and reads dataframe to calculate calculate real_n
reads <- reads %>%
  left_join(results_threshold %>% dplyr::select(X2, quantile_99_9), by = "X2")

# Subtract quantile_99_9 from n to create the new variable real_n
reads <- reads %>%
  mutate(real_n = ifelse(n - quantile_99_9 < 0, 0, n - quantile_99_9))

# Handle the specific case for X2 = intI125 where real_n should be n
reads <- reads %>%
  mutate(real_n = ifelse(X2 == "intI125" & is.na(real_n), n, real_n))

reads <- reads %>%
  mutate(real_n = ifelse(X2 == "eltB12" & is.na(real_n), n, real_n))

reads <- reads %>%
  mutate(real_n = ifelse(X2 == "eltB2318" & is.na(real_n), n, real_n))

reads[reads == 0] <- NA

#<1000
ggplot(complete(reads, X2, Sample_real), aes(X2, Sample_real)) +
  geom_tile(aes(fill = case_when(
    is.na(real_n) ~ "gray95",
    real_n > 0 & real_n <= 1 ~ "red", 
    real_n > 1 & real_n <= 100 ~ "darkorange",
    real_n > 100 & real_n <= 1000 ~ "gold1",# Assuming you want to include 1 in the range
    TRUE ~ "black"
  )), colour = "white") +
  labs(title = "Multiplex-testing on Bangladesh strains, wihtout considering PCR Neg 99.9", x = "Probe-pair", y = "Sample") +
  scale_fill_manual(values = c("gray95" = "gray95", "black" = "black", "red" = "red", "darkorange"="darkorange", "gold1"="gold1"), guide = "none") +
  theme(axis.text.x = element_text(size = 8, angle = 90), 
        axis.text.y = element_text(size = 6.5))

#Save the file-------
write.csv(reads, "bangladesh_real_n_v2.csv", row.names = F) #For final code add a step in which you make the dataframe reads complete. Otherwise, you miss false negatives, which are the most important error.
#I manually added the false positive to the file named "bangladesh_real_n.csv". 

#Re-import import reads--------
##Re-import the file "bangladesh_real_n.csv" after having added the presence/absence in WGS.
reads <- read.csv("bangladesh_real_n.csv", header = TRUE)
reads

reads[reads == 0] <- NA

# Ensure complete combinations of X2 and Sample_real
complete_reads <- complete(reads, X2, Sample_real)

# Create a new column for the fill categories
complete_reads <- complete_reads %>%
  mutate(fill_category = case_when(
    is.na(real_n) & wgs == "t" ~ "True negative",
    real_n > 0 & wgs == "y" ~ "True positive",
    real_n > 0 & wgs == "n" ~ "False positive",
    is.na(real_n) & wgs == "f" ~ "False negative",
    TRUE ~ "True negative"
  ))

# Generate the plot
ggplot(complete_reads, aes(X2, Sample_real)) +
  geom_tile(aes(fill = fill_category), colour = "white") +
  labs(title="Multiplex-testing on E. coli samples", x = "Target genes", y = "Sample") +
  scale_fill_manual(values = c(
    "True negative" = "gray95",
    "True positive" = "green",
    "False positive" = "orange",
    "False negative" = "red"
  )) +
  theme(axis.text.x = element_text(size = 8, angle = 90),
        axis.text.y = element_text(size = 5),
        legend.title = element_blank())
