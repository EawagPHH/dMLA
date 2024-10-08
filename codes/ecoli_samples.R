#Load packages--------
library (ggplot2)
library(dplyr)
library(gridExtra)
library(ggpubr)
library(tidyr)
library(seqinr)
library(readr)

#dMLA ecoli_set1 strains-------------------
##Formatting df---------
#Import results
setwd("~/switchdrive/Institution/Manuscripts/02_dMLA/dmla-amr-vfs/data")
reads <- read.csv("ecoli_set1_output.csv", header = TRUE)
reads

reads <- reads %>% 
  mutate(Sample_real = as.character(X1)) %>% 
  mutate(Sample_real = case_when(
    Sample_real == "CCACTAG" ~ "HH03C_1",
    Sample_real == "GTCTTCT" ~ "HH03CH_1",Sample_real == "ACAAAGC" ~ "HH03H_1",Sample_real == "AAAAGGC" ~ "HH03S_1",Sample_real == "ACTGTGT" ~ "HH08C_1", Sample_real == "CTGGTAC" ~ "HH08CH_1",
    Sample_real == "GTTGCTA" ~ "HH08H_1", Sample_real == "CCATTCA" ~ "zNegative_1",Sample_real == "GATATCG" ~ "HH13C_1",Sample_real == "ACAGGAT" ~ "HH13CH_1",Sample_real == "CGCATAC" ~ "HH13H_1",
    Sample_real == "GGACCTA" ~ "HH13S_1", Sample_real == "GAACTGA" ~ "HH14C_1", Sample_real == "AGTCGTG" ~ "HH14CH_1",Sample_real == "TTCCAGG" ~ "HH14H_1",Sample_real == "TAGAGCG" ~ "zNegative_2",
    Sample_real == "AAACCCT" ~ "HH14S_1",Sample_real == "ATCCAGT" ~ "HH15C_1",Sample_real == "TCTTCGT" ~ "HH15CH_1",Sample_real == "GTGTCCT" ~ "HH15H_1",Sample_real == "CACAGAT" ~ "HH16C_1",
    Sample_real == "ATAGAGC" ~ "HH16CH_1", Sample_real == "GGTCATG" ~ "HH16H_1", Sample_real == "CCATCTC" ~ "zNegative_3",Sample_real == "TATGCAG" ~ "HH16S_1",Sample_real == "AACCGCT" ~ "HH17C_1",
    Sample_real == "TGACTCT" ~ "HH17CH_1",Sample_real == "AAGCACA" ~ "HH17H_1",Sample_real == "TGTGAAC" ~ "HH17S_1",Sample_real == "AACAGTG" ~ "HH18C_1",Sample_real == "GCCTATA" ~ "HH18CH_1",
    Sample_real == "GTCCTCT" ~ "zNegative_4",Sample_real == "TAGAACG" ~ "HH18H_1",Sample_real == "TGGGTGA" ~ "HH18S_1", Sample_real == "ACACGTG" ~ "HH03C_2",Sample_real == "TGCCCAA" ~ "HH03CH_2",
    Sample_real == "TCGTACA" ~ "HH03H_2", Sample_real == "AACCAAG" ~ "HH03S_2", Sample_real == "GTGCTAT" ~ "HH08C_2",Sample_real == "TCCTCAT" ~ "zNegative_5",Sample_real == "GTCGTAT" ~ "HH08CH_2",
    Sample_real == "TCTGGAA" ~ "HH08H_2",Sample_real == "CAGTAGG" ~ "HH13C_2", Sample_real == "ACGTCAT" ~ "HH13CH_2",Sample_real == "ACGGAAC" ~ "HH13H_2",Sample_real == "CGGATGA" ~ "HH13S_2",
    Sample_real == "CCGCATA" ~ "HH14C_2", Sample_real == "GTGAAGA" ~ "zNegative_6", Sample_real == "AACGTGT" ~ "HH14CH_2",Sample_real == "AGAAGAC" ~ "HH14H_2",Sample_real == "TGTAGGG" ~ "HH14S_2",
    Sample_real == "AATGCCT" ~ "HH15C_2", Sample_real == "GCTTTCT" ~ "HH15CH_2", Sample_real == "CATGAAG" ~ "HH15H_2",Sample_real == "CGATCAC" ~ "HH16C_2",Sample_real == "TCGCGTT" ~ "zNegative_7",
    Sample_real == "GATTGGC" ~ "HH16CH_2", Sample_real == "CATCGGT" ~ "HH16H_2", Sample_real == "CTTTCCA" ~ "HH16S_2",Sample_real == "ACCAGAT" ~ "HH17C_2",Sample_real == "GTGATAC" ~ "HH17CH_2",
    Sample_real == "TGCATCC" ~ "HH17H_2", Sample_real == "TGCAAAC" ~ "HH17S_2", Sample_real == "GCAACGA" ~ "zNegative_8",Sample_real == "GCCATAC" ~ "HH18C_2",Sample_real == "TACCTTC" ~ "HH18CH_2",
    Sample_real == "GAATCGA" ~ "HH18H_2", Sample_real == "CTACGTT" ~ "HH18S_2", Sample_real == "GTTTCGG" ~ "uidA160_1",Sample_real == "GCAAATG" ~ "uidA160_2",Sample_real == "CTGACAC" ~ "zPCRNeg_1",
    Sample_real == "GCACTTT" ~ "zPCRNeg_2")) 

##Distribution of reads------------------
##Filter of negative samples containing positive
# Define the list of negative samples
negative_samples <- c("zNegative_1", "zNegative_2", "zNegative_3", "zNegative_4",
                      "zNegative_5", "zNegative_6", "zNegative_7", "zNegative_8",
                      "zPCRNeg_1", "zPCRNeg_2")

# Filter the dataframe
filtered_reads <- reads %>%
  filter(Sample_real %in% negative_samples |
           (Sample_real %in% c("uidA160_1", "uidA160_2") & X2 != "uidA160"))

#Plot distribution of false positive for each target
ggplot(filtered_reads, aes(x=n)) + 
  geom_histogram() +
  theme_minimal() +
  labs(title="Distribution of of false positive read counts for each probe-pair - Set 1 E. coli isolates",
       x="n",
       y="Frequency") +
  theme(legend.title = element_blank())+
 facet_wrap(X2 ~ ., ncol = 10, scales="free")

##Threshold calculation------------------
#Calculate the 99.9% quantile of the gamma distribution of each the false positive for each target as threshold
# Assuming 'filtered_reads' is your dataframe with columns 'n' and 'X2'
results_threshold <- filtered_reads %>%
  group_by(X2) %>%
  dplyr::summarise(
    mean_data = mean(n),
    var_data = var(n),
    .groups = 'drop'
  ) %>%
  mutate(
    shape_estimate = mean_data^2 / var_data,
    rate_estimate = mean_data / var_data,
    quantile_99_9 = qgamma(0.999, shape = shape_estimate, rate = rate_estimate)
  )

# Calculate the maximum n value for each X2 in filtered_reads
max_n_values <- filtered_reads %>%
  group_by(X2) %>%
  dplyr::summarise(max_n_in_negatives = max(n, na.rm = TRUE), .groups = 'drop')

# Replace Inf or NA values in quantile_99_9 with the highest n value from filtered_reads
results_threshold <- results_threshold %>%
  left_join(max_n_values, by = "X2") %>%
  mutate(
    quantile_99_9 = ifelse(is.infinite(quantile_99_9) | is.na(quantile_99_9), max_n_in_negatives, quantile_99_9)
  )

# Join the quantile data with the original dataframe, set quantile_99_9 to 0 for missing values
reads_with_threshold <- reads %>%
  left_join(results_threshold %>% select(X2, quantile_99_9), by = "X2") %>%
  mutate(quantile_99_9 = ifelse(is.na(quantile_99_9), 0, quantile_99_9))

# Calculate the real count
reads_with_threshold <- reads_with_threshold %>%
  mutate(real_n = n - quantile_99_9)

# Ensure real_n does not go below 0
reads_with_threshold <- reads_with_threshold %>%
  mutate(real_n = ifelse(real_n < 0, 0, real_n))

reads <- reads_with_threshold
reads[reads == 0] <- NA

##Import WGS results file-----
wgs <- read.csv("wgs_results_ecoli_set1.csv", header = TRUE)

# Remove the _1 and _2 suffixes from Sample_real in reads for comparison
reads <- reads %>%
  mutate(Sample_real_clean = sub("_[12]$", "", Sample_real))

# Create a new column in reads called 'wgs' with the specified conditions
reads <- reads %>%
  rowwise() %>%
  mutate(
    wgs = case_when(
      any(wgs$isolate == Sample_real_clean & wgs$probe == X2) & !is.na(real_n) ~ "y",   # Present in both with valid real_n
      !any(wgs$isolate == Sample_real_clean & wgs$probe == X2) & !is.na(real_n) ~ "n",  # Not present in wgs but real_n is valid
      any(wgs$isolate == Sample_real_clean & wgs$probe == X2) & is.na(real_n) ~ "f",   # Present in wgs but real_n is NA
      !any(wgs$isolate == Sample_real_clean & wgs$probe == X2) & is.na(real_n) ~ "t"                                      # Not present in wgs and real_n is NA
    )
  ) %>%
  ungroup() %>%
  select(-Sample_real_clean)

reads <- reads %>%
  mutate(
    wgs = case_when(
      Sample_real == "uidA160_1" & X2 == "uidA160" ~ "y",
      Sample_real == "uidA160_2" & X2 == "uidA160" ~ "y",
      TRUE ~ wgs  # retain existing values for other rows
    )
  )

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
set1=ggplot(complete_reads, aes(X2, Sample_real)) +
  geom_tile(aes(fill = fill_category), colour = "white") +
  labs(title="Multiplex-testing on E. coli samples - set 1", x = "Probe-pair", y = "Sample") +
  scale_fill_manual(values = c(
    "True negative" = "gray95",
    "True positive" = "palegreen3",
    "False positive" = "khaki",
    "False negative" = "lightsalmon"
  )) +
  theme(axis.text.x = element_text(size = 5, angle = 45, color = "black", hjust = 1),
        axis.text.y = element_text(size = 5,  color = "black", vjust = 0.5),
        legend.title = element_blank())

##Save the file-------
#write.csv(reads, "ecoli_set1_real_n_wgs.csv", row.names = F) 

##Count number of true/false positives/negatives---------
wgs_counts <- complete_reads %>%
  group_by(wgs) %>%
  dplyr::summarise(count = n())

# Print the counts
print(wgs_counts)

##Duplicates analyses--------
### Quantify correct probe detection among both duplicates-------
complete_reads <- complete_reads %>%
  mutate(wgs = replace_na(wgs, 't'))

complete_reads <- complete_reads %>%
  mutate(Sample_id = ifelse(grepl("^zNegative_", Sample_real), Sample_real, sub("_[12]$", "", Sample_real)))

# Group by Sample_id and X2, and then summarize the data
consistent_probes <- complete_reads %>%
  group_by(Sample_id, X2) %>%
  dplyr::summarize(consistent = n_distinct(wgs) == 1, .groups = 'drop')

# Filter to keep only those rows where the outcome is consistent
consistent_probes <- consistent_probes %>%
  filter(consistent)

# Count the number of consistent probes per Sample_id and calculate the percentage
consistent_probes_count_set1 <- consistent_probes %>%
  group_by(Sample_id) %>%
  dplyr::summarize(
    num_consistent_probes = n(),
    percent_consistent_probes = (n() / 63) * 100,  # Assuming the total number of probes is 63
    .groups = 'drop'
  )
# Print the result
print(consistent_probes_count_set1)

#write.csv(consistent_probes_count_set1, "consistent_probes_count_set1.csv", row.names = F)

#dMLA ecoli_set2 strains-------------------
##Formatting df------
#Import results
setwd("~/switchdrive/Institution/Manuscripts/02_dMLA/dmla-amr-vfs/data")
reads <- read.csv("ecoli_set2_output.csv", header = TRUE)
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

##Distribution of reads------------------
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
  labs(title="Distribution of of false positive read counts for each probe-pair - Set 2 E. coli isolates",
       x="n",
       y="Frequency") +
  theme(legend.title = element_blank())+
  facet_wrap(X2 ~ ., ncol = 10, scales="free")

##Threshold calculation------------------
#Calculate the 99.9% quantile of the gamma distribution of each the false positive for each target as threshold
# Assuming 'filtered_reads' is your dataframe with columns 'n' and 'X2'
results_threshold <- filtered_reads %>%
  group_by(X2) %>%
  dplyr::summarise(
    mean_data = mean(n),
    var_data = var(n),
    .groups = 'drop'
  ) %>%
  mutate(
    shape_estimate = mean_data^2 / var_data,
    rate_estimate = mean_data / var_data,
    quantile_99_9 = qgamma(0.999, shape = shape_estimate, rate = rate_estimate)
  )

# Calculate the maximum n value for each X2 in filtered_reads
max_n_values <- filtered_reads %>%
  group_by(X2) %>%
  dplyr::summarise(max_n_in_negatives = max(n, na.rm = TRUE), .groups = 'drop')

# Replace Inf or NA values in quantile_99_9 with the highest n value from filtered_reads
results_threshold <- results_threshold %>%
  left_join(max_n_values, by = "X2") %>%
  mutate(
    quantile_99_9 = ifelse(is.infinite(quantile_99_9) | is.na(quantile_99_9), max_n_in_negatives, quantile_99_9)
  )

# Join the quantile data with the original dataframe, set quantile_99_9 to 0 for missing values
reads_with_threshold <- reads %>%
  left_join(results_threshold %>% select(X2, quantile_99_9), by = "X2") %>%
  mutate(quantile_99_9 = ifelse(is.na(quantile_99_9), 0, quantile_99_9))

# Calculate the real count
reads_with_threshold <- reads_with_threshold %>%
  mutate(real_n = n - quantile_99_9)

# Ensure real_n does not go below 0
reads_with_threshold <- reads_with_threshold %>%
  mutate(real_n = ifelse(real_n < 0, 0, real_n))

reads <- reads_with_threshold
reads[reads == 0] <- NA

##Import WGS results file-----
wgs <- read.csv("wgs_results_ecoli_set2.csv", header = TRUE)
wgs$isolate=as.character(wgs$isolate)

# Remove the _1 and _2 suffixes from Sample_real in reads for comparison
reads <- reads %>%
  mutate(Sample_real_clean = sub("_[12]$", "", Sample_real))

# Create a new column in reads called 'wgs' with the specified conditions
reads <- reads %>%
  rowwise() %>%
  mutate(
    wgs = case_when(
      any(wgs$isolate == Sample_real_clean & wgs$probe == X2) & !is.na(real_n) ~ "y",   # Present in both with valid real_n
      !any(wgs$isolate == Sample_real_clean & wgs$probe == X2) & !is.na(real_n) ~ "n",  # Not present in wgs but real_n is valid
      any(wgs$isolate == Sample_real_clean & wgs$probe == X2) & is.na(real_n) ~ "f",   # Present in wgs but real_n is NA
      !any(wgs$isolate == Sample_real_clean & wgs$probe == X2) & is.na(real_n) ~ "t"                                      # Not present in wgs and real_n is NA
    )
  ) %>%
  ungroup() %>%
  select(-Sample_real_clean)


reads <- reads %>%
  mutate(
    wgs = case_when(
      Sample_real == "aggR11_1" & X2 == "aggR11" ~ "y",
      Sample_real == "aggR11_2" & X2 == "aggR11" ~ "y",
      Sample_real == "aaiC47_1" & X2 == "aaiC47" ~ "y",
      Sample_real == "aaiC47_2" & X2 == "aaiC47" ~ "y",
      TRUE ~ wgs  # retain existing values for other rows
    )
  )

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
set2=ggplot(complete_reads, aes(X2, Sample_real)) +
  geom_tile(aes(fill = fill_category), colour = "white") +
  labs(title="Multiplex-testing on E. coli samples - set 2", x = "Probe-pair", y = "Sample") +
  scale_fill_manual(values = c(
    "True negative" = "gray95",
    "True positive" = "palegreen3",
    "False positive" = "khaki",
    "False negative" = "lightsalmon"
  )) +
  theme(axis.text.x = element_text(size = 5, angle = 45, color = "black", hjust = 1),
        axis.text.y = element_text(size = 5,  color = "black", vjust = 0.51),
        legend.title = element_blank())

ggarrange(set1, set2, ncol=2, labels=c("A", "B"), common.legend = TRUE, legend = "bottom")

##Save the file-------
#write.csv(reads, "ecoli_set2_real_n_wgs.csv", row.names = F) 

##Count number of true/false positives/negatives---------
wgs_counts <- complete_reads %>%
  group_by(wgs) %>%
  dplyr::summarise(count = n())

# Print the counts
print(wgs_counts)

#Duplicates analyses--------
### Quantify correct probe detection among both duplicates-------
complete_reads <- complete_reads %>%
  mutate(wgs = replace_na(wgs, 't'))

complete_reads <- complete_reads %>%
  mutate(Sample_id = ifelse(grepl("^zNegative_", Sample_real), Sample_real, sub("_[12]$", "", Sample_real)))

# Group by Sample_id and X2, and then summarize the data
consistent_probes <- complete_reads %>%
  group_by(Sample_id, X2) %>%
  dplyr::summarize(consistent = n_distinct(wgs) == 1, .groups = 'drop')

# Filter to keep only those rows where the outcome is consistent
consistent_probes <- consistent_probes %>%
  filter(consistent)

# Count the number of consistent probes per Sample_id and calculate the percentage
consistent_probes_count_set2 <- consistent_probes %>%
  group_by(Sample_id) %>%
  dplyr::summarize(
    num_consistent_probes = n(),
    percent_consistent_probes = (n() / 63) * 100,  # Assuming the total number of probes is 63
    .groups = 'drop'
  )
# Print the result
print(consistent_probes_count_set2)

#write.csv(consistent_probes_count_set2, "consistent_probes_count_set2.csv", row.names = F) 

#Plot single figure------------------
#Import dataset
setwd("~/switchdrive/Institution/Manuscripts/02_dMLA/dmla-amr-vfs/data")
reads_merged <- read.csv("ecoli_merged.csv", header = TRUE)
reads_merged

# Ensure complete combinations of X2 and Sample_real
complete_reads <- complete(reads_merged, X2, Sample_real)

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
  labs(title="Multiplex-testing on E. coli samples", x = "Detected probe-pair (target gene)", y = "DNA Sample") +
  scale_fill_manual(values = c(
    "True negative" = "gray95",
    "True positive" = "palegreen3",
    "False positive" = "khaki",
    "False negative" = "lightsalmon"
  )) +
  theme(axis.text.x = element_text(size = 5, angle = 45, color = "black", hjust = 1),
        axis.text.y = element_text(size = 5,  color = "black", vjust = 0.5),
        legend.title = element_blank())

