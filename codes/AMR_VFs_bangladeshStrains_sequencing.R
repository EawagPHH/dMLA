#) Load packages
library (ggplot2)
library(dplyr)
library(gridExtra)
library(tidyr)
dev.new()

#2) Import all results (positive and negative included)
setwd("~/switchdrive/Institution/dMLA/AMRandVFsandPhylo/Eurofins/20230725_Order")
options(max.print=999999)
reads = read.table("20230725_sequencing_targets_all.txt", header= TRUE) #--->total dataset
reads

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

ggplot(reads, aes(X2,X1, fill= n))+ 
  geom_tile(colour="white")+
  scale_fill_gradient(low = "white", high = "blue")+
  theme(axis.text.x = element_text(size=5, angle=90))

ggplot(reads, aes(X2,Sample_real, fill= n))+ 
  geom_tile(colour="white")+
  scale_fill_gradient(low = "white", high = "blue")+
  theme(axis.text.x = element_text(size=5, angle=90))



##Negative controls and positive templates
#dataset - mcr_28, because I removed A1
negative = reads[reads$Sample %in% c("A1", "B2", "C4", "D6", "E8", "H12"),]

#summarize stat
s=negative %>% 
  group_by(X2) %>% 
  summarise(reads_mean = mean(n), reads_sd = sd(n))

#3) Import AGAIN all results WITH THRESHOLD(positive and negative included)
setwd("~/switchdrive/Institution/dMLA/AMRandVFsandPhylo/Eurofins/20230725_Order")
options(max.print=999999)
reads = read.table("20230725_sequencing_targets_all.txt", header= TRUE) #--->total dataset
reads

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

reads[reads == 0] <- NA

reads_filtered <- subset(reads, Sample_real != "PCR_Neg")

#Modified
#reads_filtered <- subset(reads_filtered, Sample_real != "aac6I_123")
#reads_filtered <- subset(reads_filtered, X2 != "aac6I123")

#gradient - molecular counts
ggplot(complete(reads_filtered, X2, Sample_real), aes(X2,Sample_real, fill= real_n))+ 
  geom_tile(colour="white")+
  labs(title="dMLA VFs/Phylo/ARGs Replicate n°3",
       x ="Genes Clusters", y = "Sample")+
  scale_fill_gradient(low = "lightblue", high = "blue", na.value = "gray95",  name ='Molecular counts')+
  theme(axis.text.x = element_text(size=8, angle=90), axis.text.y = element_text(size=6.5))

#unique colour - presence/absence
ggplot(complete(reads_filtered, X2, Sample_real), aes(X2, Sample_real)) +
  geom_tile(aes(fill = ifelse(is.na(real_n), "gray95", "black")), colour = "white") +
  labs(title = "dMLA VFs/Phylo/ARGs n°3", x = "Genes Clusters", y = "Sample") +
  scale_fill_manual(values = c(gray95 = "gray95", black = "black"), guide = "none") +
  theme(axis.text.x = element_text(size = 8, angle = 90), axis.text.y = element_text(size = 6.5))