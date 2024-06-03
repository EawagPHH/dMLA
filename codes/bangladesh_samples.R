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
reads <- read.csv("bangladesh_output.csv", header = TRUE)
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

#Distribution of reads------------------
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
  labs(title="Distribution of of false positive read counts for each probe-pair - Bangladesh isolates",
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

##Plot and colour-code based on the type of error
#unique colour - presence/absence
ggplot(complete(reads, X2, Sample_real), aes(X2, Sample_real)) +
  geom_tile(aes(fill = ifelse(is.na(real_n), "gray95", "black")), colour = "white") +
  labs(title="All positive", x = "Target genes", y = "Sample") +
  labs(title="Multiplex-testing on E. coli - All positive", x = "Target genes", y = "Sample") +
  scale_fill_manual(values = c(gray95 = "gray95", black = "black"), guide = "none") +
  theme(axis.text.x = element_text(size = 8, angle = 90), axis.text.y = element_text(size = 5))

# Ensure complete combinations of X2 and Sample_real
complete_reads <- complete(reads, X2, Sample_real)

# Generate the plot
ggplot(complete_reads, aes(X2, Sample_real)) +
  geom_tile(aes(fill = case_when(
    is.na(real_n) & wgs == "t" ~ "gray95",
    real_n > 0 & wgs == "y" ~ "green",
    real_n > 0 & wgs == "n" ~ "orange",
    is.na(real_n) & wgs == "f" ~ "red",
    TRUE ~ "gray95"
  )), colour = "white") +
  labs(title="Multiplex-testing on E. coli - All positive", x = "Target genes", y = "Sample") +
  scale_fill_identity(guide = "none") +
  theme(axis.text.x = element_text(size = 8, angle = 90),
        axis.text.y = element_text(size = 5))
