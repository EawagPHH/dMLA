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
setwd("~/switchdrive/Institution/Manuscripts/02_dMLA/dmla-amr-vfs/data")
reads1 <- read.csv("20240228_sequencing_targets_all_tube1.csv", header = TRUE)
reads1

reads1 <- reads1 %>% 
  mutate(Sample_real = as.character(X1)) %>% 
  mutate(Sample_real = case_when(
    Sample_real == "CCACTAG" ~ "mcr28",Sample_real == "GTCTTCT" ~ "mcr64",Sample_real == "ACAAAGC" ~ "aadA1",
    Sample_real == "AAAAGGC" ~ "aadA21",Sample_real == "ACTGTGT" ~ "aadA57", Sample_real == "CTGGTAC" ~ "aac3II63",
    Sample_real == "GTTGCTA" ~ "aac3VI23", Sample_real == "CCATTCA" ~ "zNegative1",Sample_real == "GATATCG" ~ "aac6I6",
    Sample_real == "ACAGGAT" ~ "aac6I123",Sample_real == "CGCATAC" ~ "ant2Ia1",Sample_real == "GGACCTA" ~ "aph6I29", 
    Sample_real == "GAACTGA" ~ "aph3386", Sample_real == "AGTCGTG" ~ "qnrB7",Sample_real == "TTCCAGG" ~ "qnrB35",
    Sample_real == "TAGAGCG" ~ "zNegative2",Sample_real == "AAACCCT" ~ "qnrS1",Sample_real == "ATCCAGT" ~ "ermB1",
    Sample_real == "TCTTCGT" ~ "mphA18",Sample_real == "GTGTCCT" ~ "dfrA243",Sample_real == "CACAGAT" ~ "dfrA369",
    Sample_real == "ATAGAGC" ~ "dfrA15", Sample_real == "GGTCATG" ~ "dfrA58", Sample_real == "CCATCTC" ~ "zNegative3",
    Sample_real == "TATGCAG" ~ "dfrA1",Sample_real == "AACCGCT" ~ "sul15",Sample_real == "TGACTCT" ~ "sul22",
    Sample_real == "AAGCACA" ~ "sul336",Sample_real == "TGTGAAC" ~ "tetA13",Sample_real == "AACAGTG" ~ "tetB43",
    Sample_real == "GCCTATA" ~ "tetM34",Sample_real == "GTCCTCT" ~ "zNegative4",Sample_real == "TAGAACG" ~ "intI125",
    Sample_real == "TGGGTGA" ~ "astA13", Sample_real == "ACACGTG" ~ "aatA22",Sample_real == "TGCCCAA" ~ "aatA12314",
    Sample_real == "TCGTACA" ~ "ipaH321", Sample_real == "AACCAAG" ~ "ipaH921", Sample_real == "GTGCTAT" ~ "eae111",
    Sample_real == "TCCTCAT" ~ "zNegative5",Sample_real == "GTCGTAT" ~ "eae6271",Sample_real == "TCTGGAA" ~ "stx2181",
    Sample_real == "CAGTAGG" ~ "stx23291", Sample_real == "ACGTCAT" ~ "aggR11",Sample_real == "ACGGAAC" ~ "stx1301",
    Sample_real == "CGGATGA" ~ "stx11581",Sample_real == "CCGCATA" ~ "bfpA82", Sample_real == "GTGAAGA" ~ "zNegative6",
    Sample_real == "AACGTGT" ~ "invE132",Sample_real == "AGAAGAC" ~ "eltA73",Sample_real == "TGTAGGG" ~ "eltA4510",
    Sample_real == "AATGCCT" ~ "eltB12", Sample_real == "GCTTTCT" ~ "eltB2318", Sample_real == "CATGAAG" ~ "estIa94",
    Sample_real == "CGATCAC" ~ "estIa1515",Sample_real == "TCGCGTT" ~ "zNegative7",Sample_real == "GATTGGC" ~ "bfpF56", 
    Sample_real == "CATCGGT" ~ "bfpF105", Sample_real == "CTTTCCA" ~ "aafA83",Sample_real == "ACCAGAT" ~ "aap223",
    Sample_real == "GTGATAC" ~ "aaiC17",Sample_real == "TGCATCC" ~ "aaiC47", Sample_real == "TGCAAAC" ~ "uidA160", 
    Sample_real == "GCAACGA" ~ "zNegative8",Sample_real == "GCCATAC" ~ "uidA616",Sample_real == "TACCTTC" ~ "gadph195",
    Sample_real == "GAATCGA" ~ "arpA19", Sample_real == "CTACGTT" ~ "chuA24", Sample_real == "GTTTCGG" ~ "trpA354",
    Sample_real == "GCAAATG" ~ "Tsp19",Sample_real == "CTGACAC" ~ "yjaA12", Sample_real == "GCACTTT" ~ "zPCR_neg"))

#Distribution of reads1------------------
##Filter of rows whit positive
filtered_reads1 <- reads1 %>% 
  filter(as.character(X2) != as.character(Sample_real))
#%>%  filter(X2 == "sul336")

ggplot(filtered_reads1, aes(x=n)) + 
  geom_histogram() +
  theme_minimal() +
  labs(title="Distribution of of false positive read counts for each probe-pair",
       x="n",
       y="Frequency") +
  theme(legend.title = element_blank())+
  facet_wrap(X2 ~ ., ncol = 10, scales="free")

#Calculate the 99.9% quantile of the gamma distribution of each the false positive for each target as threshold
# Assuming 'filtered_reads1' is your dataframe with columns 'n' and 'X2'
results_threshold <- filtered_reads1 %>%
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

#Join results_threshold and reads1 dataframe to calculate calculate real_n
reads1 <- reads1 %>%
  left_join(results_threshold %>% dplyr::select(X2, quantile_99_9), by = "X2")

# Subtract quantile_99_9 from n to create the new variable real_n
reads1 <- reads1 %>%
  mutate(real_n = ifelse(n - quantile_99_9 < 0, 0, n - quantile_99_9))

reads1[reads1 == 0] <- NA

#Plot results
r1=ggplot(complete(reads1, X2, Sample_real), aes(X2, Sample_real)) +
  geom_tile(aes(fill = case_when(
    is.na(real_n) ~ "n = 0",
    real_n > 0 & real_n <= 1 ~ "0 > n ≤ 1", 
    real_n > 1 & real_n <= 10 ~ "1 > n ≤ 10",
    real_n > 10 & real_n <= 100 ~ "10 > n ≤ 100",
    real_n > 100 & real_n <= 1000 ~ "100 > n ≤ 1000",# Assuming you want to include 1 in the range
    TRUE ~ "n > 1000"
  )), colour = "white") +
  labs(title = "Replicate 1", x = "Probe-pair", y = "Sample", fill = "Reads counts") +
  scale_fill_manual(values = c("n = 0" = "gray95", "0 > n ≤ 1" = "yellowgreen", "1 > n ≤ 10"="yellow2","10 > n ≤ 100"="goldenrod1", 
                               "100 > n ≤ 1000"="tomato1", "n > 1000" = "black"),
                    limits = c("n = 0", "0 > n ≤ 1", "1 > n ≤ 10", "10 > n ≤ 100", "100 > n ≤ 1000", "n > 1000")) +
  theme(axis.text.x = element_text(size = 8, angle = 90), 
        axis.text.y = element_text(size = 6.5))



#Import results tube 2
setwd("~/switchdrive/Institution/Manuscripts/02_dMLA/dmla-amr-vfs/data")
reads2 <- read.csv("20240228_sequencing_targets_all_tube2.csv", header = TRUE)
reads2

reads2 <- reads2 %>% 
  mutate(Sample_real = as.character(X1)) %>% 
  mutate(Sample_real = case_when(
    Sample_real == "CCACTAG" ~ "mcr28",Sample_real == "GTCTTCT" ~ "mcr64",Sample_real == "ACAAAGC" ~ "aadA1",
    Sample_real == "AAAAGGC" ~ "aadA21",Sample_real == "ACTGTGT" ~ "aadA57", Sample_real == "CTGGTAC" ~ "aac3II63",
    Sample_real == "GTTGCTA" ~ "aac3VI23", Sample_real == "CCATTCA" ~ "zNegative1",Sample_real == "GATATCG" ~ "aac6I6",
    Sample_real == "ACAGGAT" ~ "aac6I123",Sample_real == "CGCATAC" ~ "ant2Ia1",Sample_real == "GGACCTA" ~ "aph6I29", 
    Sample_real == "GAACTGA" ~ "aph3386", Sample_real == "AGTCGTG" ~ "qnrB7",Sample_real == "TTCCAGG" ~ "qnrB35",
    Sample_real == "TAGAGCG" ~ "zNegative2",Sample_real == "AAACCCT" ~ "qnrS1",Sample_real == "ATCCAGT" ~ "ermB1",
    Sample_real == "TCTTCGT" ~ "mphA18",Sample_real == "GTGTCCT" ~ "dfrA243",Sample_real == "CACAGAT" ~ "dfrA369",
    Sample_real == "ATAGAGC" ~ "dfrA15", Sample_real == "GGTCATG" ~ "dfrA58", Sample_real == "CCATCTC" ~ "zNegative3",
    Sample_real == "TATGCAG" ~ "dfrA1",Sample_real == "AACCGCT" ~ "sul15",Sample_real == "TGACTCT" ~ "sul22",
    Sample_real == "AAGCACA" ~ "sul336",Sample_real == "TGTGAAC" ~ "tetA13",Sample_real == "AACAGTG" ~ "tetB43",
    Sample_real == "GCCTATA" ~ "tetM34",Sample_real == "GTCCTCT" ~ "zNegative4",Sample_real == "TAGAACG" ~ "intI125",
    Sample_real == "TGGGTGA" ~ "astA13", Sample_real == "ACACGTG" ~ "aatA22",Sample_real == "TGCCCAA" ~ "aatA12314",
    Sample_real == "TCGTACA" ~ "ipaH321", Sample_real == "AACCAAG" ~ "ipaH921", Sample_real == "GTGCTAT" ~ "eae111",
    Sample_real == "TCCTCAT" ~ "zNegative5",Sample_real == "GTCGTAT" ~ "eae6271",Sample_real == "TCTGGAA" ~ "stx2181",
    Sample_real == "CAGTAGG" ~ "stx23291", Sample_real == "ACGTCAT" ~ "aggR11",Sample_real == "ACGGAAC" ~ "stx1301",
    Sample_real == "CGGATGA" ~ "stx11581",Sample_real == "CCGCATA" ~ "bfpA82", Sample_real == "GTGAAGA" ~ "zNegative6",
    Sample_real == "AACGTGT" ~ "invE132",Sample_real == "AGAAGAC" ~ "eltA73",Sample_real == "TGTAGGG" ~ "eltA4510",
    Sample_real == "AATGCCT" ~ "eltB12", Sample_real == "GCTTTCT" ~ "eltB2318", Sample_real == "CATGAAG" ~ "estIa94",
    Sample_real == "CGATCAC" ~ "estIa1515",Sample_real == "TCGCGTT" ~ "zNegative7",Sample_real == "GATTGGC" ~ "bfpF56", 
    Sample_real == "CATCGGT" ~ "bfpF105", Sample_real == "CTTTCCA" ~ "aafA83",Sample_real == "ACCAGAT" ~ "aap223",
    Sample_real == "GTGATAC" ~ "aaiC17",Sample_real == "TGCATCC" ~ "aaiC47", Sample_real == "TGCAAAC" ~ "uidA160", 
    Sample_real == "GCAACGA" ~ "zNegative8",Sample_real == "GCCATAC" ~ "uidA616",Sample_real == "TACCTTC" ~ "gadph195",
    Sample_real == "GAATCGA" ~ "arpA19", Sample_real == "CTACGTT" ~ "chuA24", Sample_real == "GTTTCGG" ~ "trpA354",
    Sample_real == "GCAAATG" ~ "Tsp19",Sample_real == "CTGACAC" ~ "yjaA12", Sample_real == "GCACTTT" ~ "zPCR_neg"))

#Distribution of reads2------------------
##Filter of rows whit positive
filtered_reads2 <- reads2 %>% 
  filter(as.character(X2) != as.character(Sample_real))
#%>%  filter(X2 == "sul336")

ggplot(filtered_reads2, aes(x=n)) + 
  geom_histogram() +
  theme_minimal() +
  labs(title="Distribution of of false positive read counts for each probe-pair",
       x="n",
       y="Frequency") +
  theme(legend.title = element_blank())+
  facet_wrap(X2 ~ ., ncol = 10, scales="free")

#Calculate the 99.9% quantile of the gamma distribution of each the false positive for each target as threshold
# Assuming 'filtered_reads2' is your dataframe with columns 'n' and 'X2'
results_threshold <- filtered_reads2 %>%
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

#Join results_threshold and reads2 dataframe to calculate calculate real_n
reads2 <- reads2 %>%
  left_join(results_threshold %>% dplyr::select(X2, quantile_99_9), by = "X2")

#Transform threshold from NaN to 1 for aaic47 and aggR11
reads2$quantile_99_9[is.na(reads2$quantile_99_9)] <- 1

# Subtract quantile_99_9 from n to create the new variable real_n
reads2 <- reads2 %>%
  mutate(real_n = ifelse(n - quantile_99_9 < 0, 0, n - quantile_99_9))

reads2[reads2 == 0] <- NA

#Plot results
r2=ggplot(complete(reads2, X2, Sample_real), aes(X2, Sample_real)) +
  geom_tile(aes(fill = case_when(
    is.na(real_n) ~ "n = 0",
    real_n > 0 & real_n <= 1 ~ "0 > n ≤ 1", 
    real_n > 1 & real_n <= 10 ~ "1 > n ≤ 10",
    real_n > 10 & real_n <= 100 ~ "10 > n ≤ 100",
    real_n > 100 & real_n <= 1000 ~ "100 > n ≤ 1000",# Assuming you want to include 1 in the range
    TRUE ~ "n > 1000"
  )), colour = "white") +
  labs(title = "Replicate 2", x = "Probe-pair", y = "Sample", fill = "Reads counts") +
  scale_fill_manual(values = c("n = 0" = "gray95", "0 > n ≤ 1" = "yellowgreen", "1 > n ≤ 10"="yellow2","10 > n ≤ 100"="goldenrod1", 
                               "100 > n ≤ 1000"="tomato1", "n > 1000" = "black"),
                    limits = c("n = 0", "0 > n ≤ 1", "1 > n ≤ 10", "10 > n ≤ 100", "100 > n ≤ 1000", "n > 1000")) +
  theme(axis.text.x = element_text(size = 8, angle = 90), 
        axis.text.y = element_text(size = 6.5))

ggarrange(r1, r2, ncol=2, labels=c("A", "B"), common.legend = TRUE, legend = "right")
