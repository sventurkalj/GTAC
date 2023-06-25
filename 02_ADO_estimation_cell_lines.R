#Estimation of ADO in cell lines  
#If possible, it is a good idea to estimate allelic drop-out (ADO) for an amplicon of interest, if there are cell lines or models containing a heterozygous variant at that locus

#Load relevant libraries
library(ggpubr)
library(tidyverse)
library(stringr)
library(cowplot)
theme_set(theme_cowplot())

#Genotype calling - this part of the code is exactly the same one used in the previous section
#### Read in files generated with the pre-processing pipeline; a separate file for each genotyped locus is created ####
genotyping_gDNA <- read_tsv('locus_n.txt')

#### Subtract the file for cells of interest and name the samples ####
genotyping_gDNA <- genotyping_gDNA[c(1:192),] %>% # Get only relevant cells for locus
  mutate(Sample = case_when(str_detect(id, 'K562_Tn5') ~ 'K562',
                            str_detect(id, 'THP1_Tn5') ~ 'THP1',
                            TRUE ~ 'NA')
  )     

# Name of control sample
control <- "THP1"

# Sample being analysed
sample <- "K562"

#### Define mutation and WT call ####
# This sets which columns in the file to use as reference and variant alleles
ref <- quo(T)
mut <- quo(G)

#### Calculate coverage and genotyping success rate####
gDNA_min_coverage <- 30

genotyping_gDNA <- genotyping_gDNA %>%
  mutate(coverage = A+C+G+T+del+ins) %>%
  mutate(captured = case_when(coverage >= gDNA_min_coverage ~ 'yes',
                              TRUE ~ 'no'))

# Understand how many cells show efficient amplification
success <- genotyping_gDNA %>% 
  group_by(Sample, captured) %>%
  dplyr::count() 

### Filter the data set to remove failed libraries and to filter for samples of interest ###
genotyping_gDNA_filtered <- genotyping_gDNA %>% 
  filter(coverage >= gDNA_min_coverage & Sample %in% c(sample, control)) 

###### Calculate scVAF and call genotypes ####

genotyping_gDNA_filtered <- genotyping_gDNA_filtered %>%
  mutate(VAF=!!mut/coverage)

# Get control values
genotyping_gDNA_filtered %>% 
  filter(Sample == control) %>% 
  summarise(mean_VAF = mean(VAF),
            max_VAF = max(VAF),
            min_VAF = min(VAF),
            mean_Mut_reads = mean(!!mut),
            max_Mut_reads = max(!!mut),
            SD = sd(VAF),
            '2XSD' = 2*SD,
            'Max+2XSD' = max_VAF + (2*SD)
  )

# Set the VAF threshold for calling a mutant cell based on these values

gDNA_min_VAF <- 0.00299 # This is max+VAF + 2XSD
gDNA_max_VAF <- 1 - gDNA_min_VAF # 1 - gDNA_min_VAF -> it is the inverse threshold
gDNA_min_mut_reads <- 9 # Max number of mutant reads in control

# Call genotype

genotyping_gDNA_filtered <- genotyping_gDNA_filtered %>% 
  mutate(genotype = case_when(
    VAF > gDNA_min_VAF & !!mut >= gDNA_min_mut_reads ~ "Mut",
    TRUE ~ "WT")) %>%
  write_csv("locus_n_integrate.csv")

#### ADO estimate - how many K562 cells are called homozygous WT or homozygous mutant when they should be heterozygous mutant? ####

genotyping_gDNA_filtered <- genotyping_gDNA_filtered %>%
  filter(Sample == 'K562') %>%
  mutate(WT_ADO = case_when(
    VAF >= gDNA_max_VAF ~ 'Yes',
    TRUE ~ 'No'
  )) %>%
  mutate(Mut_ADO = case_when(
    VAF <= gDNA_min_VAF | !!mut < gDNA_min_mut_reads ~ 'Yes',
    TRUE ~ 'No'
  )) %>%
  mutate(biallelic_detection = case_when(
    WT_ADO == 'No' & Mut_ADO == 'No' ~ 'Yes',
    TRUE ~ 'No'
  ))

total_genotyped <- nrow(genotyping_gDNA_filtered)

ado <- genotyping_gDNA_filtered %>% 
  summarise(mut_ADO = sum(Mut_ADO == 'Yes'),
            mut_ADO_percentage = (mut_ADO/total_genotyped)*100,
            wt_ADO = sum(WT_ADO == 'Yes'),
            wt_ADO_percentage = (wt_ADO/total_genotyped)*100,
            allelic_bias = sum(biallelic_detection == 'No'),
            overall_bias_percentage = (allelic_bias/total_genotyped)*100,
            average_VAF = mean(VAF)
  ) %>%
  write_csv('K562_locus_n_ADO_values.csv')

# This part identifies which K562 cells have achieved mutant detection for all alleles
# First merge genotype call files for all loci for all cells - all_mutations_all_cells.csv

k562 <- read_csv('all_mutations_all_cells.csv') %>%
  filter(Sample == "K562") %>%
  mutate(biallelic_all = case_when(
    ASXL1 == 'Mut' & HOXA9 == 'Mut' & NOTCH1 == 'Mut' & CBL == 'Mut' & TET2 == 'Mut' ~ 'Yes',
    TRUE ~ 'No'
  ))

k562 <- k562 %>%
  mutate(biallelic_4 = case_when(
    ASXL1 == 'Mut' & HOXA9 == 'Mut' & NOTCH1 == 'Mut' & CBL == 'Mut' & TET2 == 'WT' ~ 'Yes',
    ASXL1 == 'Mut' & HOXA9 == 'Mut' & NOTCH1 == 'Mut' & CBL == 'WT' & TET2 == 'Mut' ~ 'Yes',
    ASXL1 == 'Mut' & HOXA9 == 'Mut' & NOTCH1 == 'WT' & CBL == 'Mut' & TET2 == 'Mut' ~ 'Yes',
    ASXL1 == 'Mut' & HOXA9 == 'WT' & NOTCH1 == 'Mut' & CBL == 'Mut' & TET2 == 'Mut' ~ 'Yes',
    ASXL1 == 'WT' & HOXA9 == 'Mut' & NOTCH1 == 'Mut' & CBL == 'Mut' & TET2 == 'Mut' ~ 'Yes',
    TRUE ~ 'No'
  )) %>%
  write_csv("K562_geno_QC.csv")

# Get genotyping QC values based on multiplexed detection of mutations (Total of 240 cells screened)

qc <- k562 %>% 
  summarise(all_detected = sum(biallelic_all == 'Yes'),
            percentage_all_detected = all_detected/240,
            four_detected = sum(biallelic_4 == 'Yes'),
            percentage_three_detected = three_detected/240
  ) %>%
  mutate(less_than_four_detected = 1 - (percentage_all_detected + percentage_three_detected))

# Understand number of cells with successful genotype detection and successful scATAC QC
# scATAC-seq pre-processing pipeline will give as output the sample_info.csv file

atac <- read_csv('sample_info.csv')
geno <- read_csv('K562_geno_QC.csv')
atac_geno <- left_join(atac,geno,by=cell_name)

success_rate <- atac_geno %>%
  mutate(initial_qc = case_when(
    mapping_rate >= 90 & uniq_nuc_frags >= 10000 &  biallelic_all == 'yes' ~ 'Passed_and_all_loci',
    mapping_rate >= 90 & uniq_nuc_frags >= 10000 &  biallelic_4 == 'yes' ~ 'Passed_and_four_loci',
    TRUE ~ 'Failed'
  ))
