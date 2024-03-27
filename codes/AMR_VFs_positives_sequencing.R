#dMLA Positive Templates-------------------
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
reads <- read.csv("20230725_sequencing_targets_all_noPCRneg.csv", header = TRUE)
reads

#Assign barcodes to samples
reads <- reads %>% 
  mutate(Sample = as.character(X1)) %>% 
  mutate(Sample = case_when(
    Sample == "GAACACA" ~ "A1", Sample == "ACATTCA" ~ "A2",Sample == "TCGCTAG" ~ "A3", Sample == "ATCATTA" ~ "A4",Sample == "TGTATGT" ~ "A5",
    Sample == "TAAGATA" ~ "A6",Sample == "CGTTTCA" ~ "A7",Sample == "ACGTTGC" ~ "A8",Sample == "TAGATGA" ~ "A9", Sample == "CCACTAG" ~ "A10",
    Sample == "GTCTTCT" ~ "A11",Sample == "ACAAAGC" ~ "A12",Sample == "AAAAGGC" ~ "B1",Sample == "ACTGTGT" ~ "B2", Sample == "CTGGTAC" ~ "B3",
    Sample == "GTTGCTA" ~ "B4", Sample == "CCATTCA" ~ "B5",Sample == "GATATCG" ~ "B6",Sample == "ACAGGAT" ~ "B7",Sample == "CGCATAC" ~ "B8",
    Sample == "GGACCTA" ~ "B9", Sample == "GAACTGA" ~ "B10", Sample == "AGTCGTG" ~ "B11",Sample == "TTCCAGG" ~ "B12",Sample == "TAGAGCG" ~ "C1",
    Sample == "AAACCCT" ~ "C2",Sample == "ATCCAGT" ~ "C3",Sample == "TCTTCGT" ~ "C4",Sample == "GTGTCCT" ~ "C5",Sample == "CACAGAT" ~ "C6",
    Sample == "ATAGAGC" ~ "C7", Sample == "GGTCATG" ~ "C8", Sample == "CCATCTC" ~ "C9",Sample == "TATGCAG" ~ "C10",Sample == "AACCGCT" ~ "C11",
    Sample == "TGACTCT" ~ "C12",Sample == "AAGCACA" ~ "D1",Sample == "TGTGAAC" ~ "D2",Sample == "AACAGTG" ~ "D3",Sample == "GCCTATA" ~ "D4",
    Sample == "GTCCTCT" ~ "D5",Sample == "TAGAACG" ~ "D6",Sample == "TGGGTGA" ~ "D7", Sample == "ACACGTG" ~ "D8",Sample == "TGCCCAA" ~ "D9",
    Sample == "TCGTACA" ~ "D10", Sample == "AACCAAG" ~ "D11", Sample == "GTGCTAT" ~ "D12",Sample == "TCCTCAT" ~ "E1",Sample == "GTCGTAT" ~ "E2",
    Sample == "TCTGGAA" ~ "E3",Sample == "CAGTAGG" ~ "E4", Sample == "ACGTCAT" ~ "E5",Sample == "ACGGAAC" ~ "E6",Sample == "CGGATGA" ~ "E7",
    Sample == "CCGCATA" ~ "E8", Sample == "GTGAAGA" ~ "E9", Sample == "AACGTGT" ~ "E10",Sample == "AGAAGAC" ~ "E11",Sample == "TGTAGGG" ~ "E12",
    Sample == "AATGCCT" ~ "F1", Sample == "GCTTTCT" ~ "F2", Sample == "CATGAAG" ~ "F3",Sample == "CGATCAC" ~ "F4",Sample == "TCGCGTT" ~ "F5",
    Sample == "GATTGGC" ~ "F6", Sample == "CATCGGT" ~ "F7", Sample == "CTTTCCA" ~ "F8",Sample == "ACCAGAT" ~ "F9",Sample == "GTGATAC" ~ "F10",
    Sample == "TGCATCC" ~ "F11", Sample == "TGCAAAC" ~ "F12", Sample == "GCAACGA" ~ "G1",Sample == "GCCATAC" ~ "G2",Sample == "TACCTTC" ~ "G3",
    Sample == "GAATCGA" ~ "G4", Sample == "CTACGTT" ~ "G5", Sample == "GTTTCGG" ~ "G6",Sample == "GCAAATG" ~ "G7",Sample == "CTGACAC" ~ "G8",
    Sample == "GCACTTT" ~ "G9", Sample == "AGCAAAC" ~ "G10", Sample == "ATCCTCG" ~ "G11",Sample == "CTGCTGT" ~ "G12",Sample == "CAATCAC" ~ "H1",
    Sample == "GTTCTGG" ~ "H2", Sample == "TTGCTAC" ~ "H3", Sample == "CATACGT" ~ "H4",Sample == "GATCGTT" ~ "H5",Sample == "CACGTTT" ~ "H6",
    Sample == "AACACCT" ~ "H7", Sample == "AGCTGTA" ~ "H8", Sample == "GGTCTTG" ~ "H9",Sample == "CCGACAT" ~ "H10",Sample == "TATCGTC" ~ "H11",
    Sample == "ATGGTTG" ~ "H12"))

reads <- reads %>% 
  mutate(Sample_real = as.character(X1)) %>% 
  mutate(Sample_real = case_when(
    Sample_real == "GAACACA" ~ "empty", Sample_real == "ACATTCA" ~ "HH03C_1",Sample_real == "TCGCTAG" ~ "HH03CH_1", Sample_real == "ATCATTA" ~ "HH03H_1",Sample_real == "TGTATGT" ~ "HH03S_1",
    Sample_real == "TAAGATA" ~ "HH08C_1",Sample_real == "CGTTTCA" ~ "HH08CH_1",Sample_real == "ACGTTGC" ~ "HH08H_1",Sample_real == "TAGATGA" ~ "HH13C_1", Sample_real == "CCACTAG" ~ "HH13CH_1",
    Sample_real == "GTCTTCT" ~ "HH13H_1",Sample_real == "ACAAAGC" ~ "HH13S_1",Sample_real == "AAAAGGC" ~ "HH14C_1",Sample_real == "ACTGTGT" ~ "Negative1", Sample_real == "CTGGTAC" ~ "HH14CH_1",
    Sample_real == "GTTGCTA" ~ "HH14H_1", Sample_real == "CCATTCA" ~ "HH14S_1",Sample_real == "GATATCG" ~ "HH15C_1",Sample_real == "ACAGGAT" ~ "HH15CH_1",Sample_real == "CGCATAC" ~ "HH15H_1",
    Sample_real == "GGACCTA" ~ "HH16C_1", Sample_real == "GAACTGA" ~ "HH16CH_1", Sample_real == "AGTCGTG" ~ "HH16H_1",Sample_real == "TTCCAGG" ~ "HH16S_1",Sample_real == "TAGAGCG" ~ "HH17C_1",
    Sample_real == "AAACCCT" ~ "HH17CH_1",Sample_real == "ATCCAGT" ~ "HH17H_1",Sample_real == "TCTTCGT" ~ "Negative2",Sample_real == "GTGTCCT" ~ "HH17S_1",Sample_real == "CACAGAT" ~ "HH18C_1",
    Sample_real == "ATAGAGC" ~ "HH18CH_1", Sample_real == "GGTCATG" ~ "HH18H_1", Sample_real == "CCATCTC" ~ "HH18S_1",Sample_real == "TATGCAG" ~ "HH03C_2",Sample_real == "AACCGCT" ~ "HH03CH_2",
    Sample_real == "TGACTCT" ~ "HH03H_2",Sample_real == "AAGCACA" ~ "HH03S_2",Sample_real == "TGTGAAC" ~ "HH08C_2",Sample_real == "AACAGTG" ~ "HH08CH_2",Sample_real == "GCCTATA" ~ "HH08H_2",
    Sample_real == "GTCCTCT" ~ "HH13C_2",Sample_real == "TAGAACG" ~ "Negative3",Sample_real == "TGGGTGA" ~ "HH13CH_2", Sample_real == "ACACGTG" ~ "HH13H_2",Sample_real == "TGCCCAA" ~ "HH13S_2",
    Sample_real == "TCGTACA" ~ "HH14C_2", Sample_real == "AACCAAG" ~ "HH14CH_2", Sample_real == "GTGCTAT" ~ "HH14H_2",Sample_real == "TCCTCAT" ~ "HH14S_2",Sample_real == "GTCGTAT" ~ "HH15C_2",
    Sample_real == "TCTGGAA" ~ "HH15CH_2",Sample_real == "CAGTAGG" ~ "HH15H_2", Sample_real == "ACGTCAT" ~ "HH16C_2",Sample_real == "ACGGAAC" ~ "HH16CH_2",Sample_real == "CGGATGA" ~ "HH16H_2",
    Sample_real == "CCGCATA" ~ "Negative4", Sample_real == "GTGAAGA" ~ "HH16S_2", Sample_real == "AACGTGT" ~ "HH17C_2",Sample_real == "AGAAGAC" ~ "HH17CH_2",Sample_real == "TGTAGGG" ~ "HH17H_2",
    Sample_real == "AATGCCT" ~ "HH17S_2", Sample_real == "GCTTTCT" ~ "HH18C_2", Sample_real == "CATGAAG" ~ "HH18CH_2",Sample_real == "CGATCAC" ~ "HH18H_2",Sample_real == "TCGCGTT" ~ "HH18S_2",
    Sample_real == "GATTGGC" ~ "HH03C_3", Sample_real == "CATCGGT" ~ "HH03CH_3", Sample_real == "CTTTCCA" ~ "HH03H_3",Sample_real == "ACCAGAT" ~ "HH03S_3",Sample_real == "GTGATAC" ~ "HH08C_3",
    Sample_real == "TGCATCC" ~ "HH08CH_3", Sample_real == "TGCAAAC" ~ "HH08H_3", Sample_real == "GCAACGA" ~ "HH13C_3",Sample_real == "GCCATAC" ~ "HH13CH_3",Sample_real == "TACCTTC" ~ "HH13H_3",
    Sample_real == "GAATCGA" ~ "HH13S_3", Sample_real == "CTACGTT" ~ "HH14C_3", Sample_real == "GTTTCGG" ~ "HH14CH_3",Sample_real == "GCAAATG" ~ "HH14H_3",Sample_real == "CTGACAC" ~ "HH14S_3",
    Sample_real == "GCACTTT" ~ "HH15C_3", Sample_real == "AGCAAAC" ~ "HH15CH_3", Sample_real == "ATCCTCG" ~ "HH15H_3",Sample_real == "CTGCTGT" ~ "HH16C_3",Sample_real == "CAATCAC" ~ "HH16CH_3",
    Sample_real == "GTTCTGG" ~ "HH16H_3", Sample_real == "TTGCTAC" ~ "HH16S_3", Sample_real == "CATACGT" ~ "HH17C_3",Sample_real == "GATCGTT" ~ "HH17CH_3",Sample_real == "CACGTTT" ~ "HH17H_3",
    Sample_real == "AACACCT" ~ "HH17S_3", Sample_real == "AGCTGTA" ~ "HH18C_3", Sample_real == "GGTCTTG" ~ "HH18CH_3",Sample_real == "CCGACAT" ~ "HH18H_3",Sample_real == "TATCGTC" ~ "HH18S_3",
    Sample_real == "ATGGTTG" ~ "PCR_Neg"))

##Calculate detection limit-----------------------
# Create a reference dataframe with all levels of X2 present in the reads df
reads <- subset(reads, Sample_real != "PCR_Neg") #Remove this samples because contamination occurred and I'm aware of it.
complete_X2 <- data.frame(X2 = levels(factor(reads$X2)))

#Create a dataset with number of reads in negative controls of the reaction
neg = reads[reads$Sample_real %in% c("Negative1", "Negative2", "Negative3", "Negative4"), ]

#Compute the detection limit of the reads in the negative through the sum of the mean + 3 times the sd.
summary_stat = neg %>% 
  group_by(X2) %>% 
  summarise(reads_mean = mean(n), reads_sd = sd(n))%>%
  mutate(threshold = reads_mean + 3 * reads_sd)

# Perform left join to include all levels in summary_stat and fill other columns with NA and threshold=0
summary_stat <- left_join(complete_X2, summary_stat, by = "X2") %>%
  mutate(reads_mean = ifelse(is.na(reads_mean), NA, reads_mean),
         reads_sd = ifelse(is.na(reads_sd), NA, reads_sd),
         threshold = ifelse(is.na(threshold), NA, threshold)) %>%
  replace_na(list(threshold = 0))

# Left join reads with summary_stat on X2
reads <- left_join(reads, summary_stat, by = "X2")

# Calculate real_n based on the formula n - threshold, ensuring it is not lower than 0 and remove variables reads_mean and reads_sd
reads <- reads %>%
  mutate(real_n = ifelse(is.na(threshold), n, pmax(n - threshold, 0)))%>%
  select(-reads_mean, -reads_sd)

##Plot results------------
reads[reads == 0] <- NA

#gradient - molecular counts
ggplot(complete(reads, X2, Sample_real), aes(X2,Sample_real, fill= real_n))+ 
  geom_tile(colour="white")+
  labs(title="dMLA VFs/Phylo/ARGs Replicate n°3",
       x ="Genes Clusters", y = "Sample")+
  scale_fill_gradient(low = "lightblue", high = "blue", na.value = "gray95",  name ='Molecular counts')+
  theme(axis.text.x = element_text(size=8, angle=90), axis.text.y = element_text(size=6.5))

#unique colour - presence/absence
all=ggplot(complete(reads, X2, Sample_real), aes(X2, Sample_real)) +
  geom_tile(aes(fill = ifelse(is.na(real_n), "gray95", "black")), colour = "white") +
  labs(title="Multiplex-testing on E. coli - All positive", x = "Target genes", y = "Sample") +
  scale_fill_manual(values = c(gray95 = "gray95", black = "black"), guide = "none") +
  theme(axis.text.x = element_text(size = 8, angle = 90), axis.text.y = element_text(size = 5))

#< 1
ggplot(complete(reads, X2, Sample_real), aes(X2, Sample_real)) +
  geom_tile(aes(fill = case_when(
    is.na(real_n) ~ "gray95",
    real_n > 0 & real_n <= 1 ~ "red", # Assuming you want to include 1 in the range
    TRUE ~ "black"
  )), colour = "white") +
  labs(title = "Multiplex-testing - Tube 2", x = "Probe-pair", y = "Sample") +
  scale_fill_manual(values = c("gray95" = "gray95", "black" = "black", "red" = "red"), guide = "none") +
  theme(axis.text.x = element_text(size = 8, angle = 90), 
        axis.text.y = element_text(size = 6.5))

#<100
ggplot(complete(reads, X2, Sample_real), aes(X2, Sample_real)) +
  geom_tile(aes(fill = case_when(
    is.na(real_n) ~ "gray95",
    real_n > 0 & real_n <= 1 ~ "red", 
    real_n > 1 & real_n <= 100 ~ "darkorange",  # Assuming you want to include 1 in the range
    TRUE ~ "black"
  )), colour = "white") +
  labs(title = "Multiplex-testing - Tube 2", x = "Probe-pair", y = "Sample") +
  scale_fill_manual(values = c("gray95" = "gray95", "black" = "black", "red" = "red", "darkorange"="darkorange"), guide = "none") +
  theme(axis.text.x = element_text(size = 8, angle = 90), 
        axis.text.y = element_text(size = 6.5))

#<1000
ggplot(complete(reads, X2, Sample_real), aes(X2, Sample_real)) +
  geom_tile(aes(fill = case_when(
    is.na(real_n) ~ "gray95",
    real_n > 0 & real_n <= 1 ~ "red", 
    real_n > 1 & real_n <= 100 ~ "darkorange",
    real_n > 100 & real_n <= 1000 ~ "gold1",# Assuming you want to include 1 in the range
    TRUE ~ "black"
  )), colour = "white") +
  labs(title = "Multiplex-testing - Tube 2", x = "Probe-pair", y = "Sample") +
  scale_fill_manual(values = c("gray95" = "gray95", "black" = "black", "red" = "red", "darkorange"="darkorange", "gold1"="gold1"), guide = "none") +
  theme(axis.text.x = element_text(size = 8, angle = 90), 
        axis.text.y = element_text(size = 6.5))

## Export data--------
write.csv(reads, file = "reads_counts.csv", row.names = FALSE)

## Replicates analyses-----
# Remove negative controls
reads <- subset(reads, !(Sample_real %in% c("Negative1", "Negative2", "Negative3", "Negative4")))

# Remove rows where real_n is NA
reads <- subset(reads, !is.na(real_n))

# Extract the common part of Sample_real names before the "_"
reads <- reads %>%
  mutate(Sample_real = sub("_.*", "", Sample_real))

# Count the number of replicates for each Sample_name within each X2 group
replicate_counts <- reads %>%
  group_by(X2, Sample_real) %>%
  summarise(replicate_count = n()) %>%
  ungroup()

### Targets present in all three replicates---------
# Filter the data frame to keep only rows where all replicates are present for a unique X2
tri_reads <- reads %>%
  inner_join(replicate_counts, by = c("X2", "Sample_real")) %>%
  group_by(X2) %>%
  filter((replicate_count >= 3)) %>%
  ungroup() %>%
  dplyr::select(X2, Sample_real, real_n)

#Plot presence/absence
tri=ggplot(complete(tri_reads, X2, Sample_real), aes(X2, Sample_real)) +
  geom_tile(aes(fill = ifelse(is.na(real_n), "gray95", "black")), colour = "white") +
  labs(title = "Multiplex-testing on E. coli - Triplicates positive only", x = "Target genes", y = "Sample") +
  scale_fill_manual(values = c(gray95 = "gray95", black = "black"), guide = "none") +
  theme(axis.text.x = element_text(size = 8, angle = 90), axis.text.y = element_text(size = 5))

#<1000
ggplot(complete(tri_reads, X2, Sample_real), aes(X2, Sample_real)) +
  geom_tile(aes(fill = case_when(
    is.na(real_n) ~ "gray95",
    real_n > 0 & real_n <= 1 ~ "red", 
    real_n > 1 & real_n <= 100 ~ "darkorange",
    real_n > 100 & real_n <= 1000 ~ "gold1",# Assuming you want to include 1 in the range
    TRUE ~ "black"
  )), colour = "white") +
  labs(title = "Multiplex-testing on E. coli - Triplicates positive only", x = "Probe-pair", y = "Sample") +
  scale_fill_manual(values = c("gray95" = "gray95", "black" = "black", "red" = "red", "darkorange"="darkorange", "gold1"="gold1"), guide = "none") +
  theme(axis.text.x = element_text(size = 8, angle = 90), 
        axis.text.y = element_text(size = 5))

### Targets present in two replicates---------
# Filter the data frame to keep only rows where all replicates are present for a unique X2
dup_reads <- reads %>%
  inner_join(replicate_counts, by = c("X2", "Sample_real")) %>%
  group_by(X2) %>%
  filter((replicate_count >= 2)) %>%
  ungroup() %>%
  dplyr::select(X2, Sample_real, real_n)

#Plot presence/absence
dup=ggplot(complete(dup_reads, X2, Sample_real), aes(X2, Sample_real)) +
  geom_tile(aes(fill = ifelse(is.na(real_n), "gray95", "black")), colour = "white") +
  labs(title = "Multiplex-testing on E. coli - Duplicates positive only", x = "Target genes", y = "Sample") +
  scale_fill_manual(values = c(gray95 = "gray95", black = "black"), guide = "none") +
  theme(axis.text.x = element_text(size = 8, angle = 90), axis.text.y = element_text(size = 5))

#Distribution of reads------------------
##Filter of rows whit positive
filtered_reads = reads[reads$Sample_real %in% c("Negative1", "Negative2", "Negative3", "Negative4","PCR_Neg"), ]

#Plot distribution of false positive for each target
ggplot(filtered_reads, aes(x=n)) + 
  geom_histogram() +
  theme_minimal() +
  labs(title="Distribution of of false positive read counts for each probe-pair (Tube2)",
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
