# Analysis of GTAC single cell genotyping data

#### Load the relevant libraries ####
library(ggpubr)
library(tidyverse)
library(stringr)
library(cowplot)
theme_set(theme_cowplot())

#### Read in files generated with the pre-processing pipeline; a separate file for each genotyped locus is created ####
genotyping_gDNA <- read_tsv('locus_n.txt')

#### Subtract the file for cells of interest and name the samples ####
genotyping_gDNA <- genotyping_gDNA[c(1:192),] %>% # Get only relevant cells for locus
  mutate(Sample = case_when(str_detect(id, 'AML_Tn5') ~ 'AML',
                            str_detect(id, 'WT_BM_Tn5') ~ 'WT_BM',
                            TRUE ~ 'NA')
  )     

# Name of control sample
control <- "WT_BM"

# Sample being analysed
sample <- "AML"

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
gDNA_max_VAF <- 1 - gDNA_min_VAF # 1 - gDNA_min_VAD -> it is the inverse threshold
gDNA_min_mut_reads <- 9 # Max number of mutant reads in control

# Call genotype

genotyping_gDNA_filtered <- genotyping_gDNA_filtered %>% 
  mutate(genotype = case_when(
    VAF > gDNA_min_VAF & !!mut >= gDNA_min_mut_reads ~ "Mut",
    TRUE ~ "WT")) %>%
  write_csv("locus_n_integrate.csv")

# Plot mutant reads vs scVAF - thresholds for ADO estimation

ggplot(genotyping_gDNA_filtered, aes(x = !!mut, y = VAF, color = Sample)) +
  geom_rect(aes(xmin = gDNA_min_mut_reads, xmax = Inf,  ymin = gDNA_min_VAF, ymax = gDNA_max_VAF),
            fill = "azure", colour = NA) + geom_point(position = "jitter")+
  geom_point(size=1.5) +
  geom_vline(aes(xintercept = gDNA_min_mut_reads),  linetype = "longdash")+
  geom_hline(aes(yintercept = gDNA_min_VAF), linetype = "longdash") +
  geom_hline(aes(yintercept = gDNA_max_VAF), linetype = "longdash", color = 'red') +
  labs(title = 'Locus n', y = "scVAF", x = "Mutant reads") +
  scale_x_log10()
ggsave("Thresholds_locus_n.pdf", device = "pdf", width = 7)

###### Repeat this for each amplified locus before proceeding to genotype integration ####
