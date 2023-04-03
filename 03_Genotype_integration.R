# Integrate genotypes from GTAC single cell genotyping data

#### Load the relevant libraries ####
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
library(forcats)
library(readr)
library(dplyr)

# Name of control sample
control <- "control_BM"

# Sample being analysed
sample <- "AML"

##### READ IN FILES - these files are produced after genotypes are called for each locus #####

mutation_1_path <- "TET2_integrate.csv"
mutation_2_path <- "ASXL1_integrate.csv"
mutation_3_path <- "RUNX1_integrate.csv"
mutation_4_path <- "BCOR_1_integrate.csv"
mutation_5_path <- "BCOR_2_integrate.csv"

# Define a function that reads file and takes columns needed

read_genotyping_file <- function(file_path) {
  file <- read_csv(file_path) %>% 
    select(cell_id, genotype, mutation)
  mutation_name <- file$mutation[1]
  file <- file %>% rename(!!mutation_name := genotype) %>% 
    select(-mutation)
  return(file)
}

# Read in the genotyping data

mutation_1 <- read_genotyping_file(mutation_1_path)
mutation_2 <- read_genotyping_file(mutation_2_path)
mutation_3 <- read_genotyping_file(mutation_3_path)
mutation_4 <- read_genotyping_file(mutation_4_path)
mutation_5 <- read_genotyping_file(mutation_5_path)

genotype_list = list(mutation_1, mutation_2, mutation_3, mutation_4, mutation_5)
genotypes <- genotype_list %>% reduce(inner_join, by='cell_id')

##### Generate a raster plot to perform visual inspection of the clonal hierarchy first #####

genotypes <- genotypes %>% 
  mutate(ASXL1_13_mut = factor(ASXL1_13_mut, levels = c("Mut", "WT", "Undetected")),
         TET2_6_mut1 = factor(TET2_6_mut1, levels = c("Mut", "WT","Undetected")),
         RUX1_5_mut = factor(RUX1_5_mut, levels = c("Mut", "WT","Undetected")),
         BCOR_4_insertion1 = factor(BCOR_4_insertion1, levels = c("Mut", "WT","Undetected")),
         BCOR_4_insertion2_1 = factor(BCOR_4_insertion2_1, levels = c("Mut", "WT","Undetected"))) %>%
  arrange(TET2_6_mut1, ASXL1_13_mut, RUX1_5_mut, BCOR_4_insertion1, BCOR_4_insertion2_1) %>% 
  mutate(cell_id=factor(cell_id, levels=unique(cell_id)))

genotypes_tidied <- genotypes %>% 
  gather(TET2_6_mut1, ASXL1_13_mut, RUX1_5_mut, BCOR_4_insertion1, BCOR_4_insertion2_1, key = "Mutation", value = "genotype", factor_key = T) %>% 
  mutate(genotype = factor(genotype, levels = c("Mut", "WT", "Undetected"))) 

p1 <- ggplot(genotypes_tidied,
             aes(x =  cell_id, y = Mutation, fill=genotype))+
  geom_tile()+
  geom_vline(xintercept=genotypes_tidied$cell_id, size=0.005, colour="grey50")+
  scale_fill_manual(values = c("WT" = "#bdbdbd", "Mut" = "cyan3", "Undetected" = "#f0f0f0"))+ 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line = element_blank(),
        axis.text.y = element_text(size = 16), axis.title.y = element_text(size = 16))+
  labs(title = sample, x = NULL)
p1

##### Assign clonal identity; This is based on the most likely identity from the raster plot. Note that this structure has to be validated with infSCITE as stated in STAR Methods #####
# This part of the code needs to be adapted for every individual clonal hierarchy - the clonal hierarchy here is the one presented in the paper

genotypes <- genotypes %>% 
  mutate(clone = case_when(
    ASXL1_13_mut == 'Mut' & BCOR_4_insertion1 == 'Mut' & BCOR_4_insertion2_1 == 'Mut' & RUX1_5_mut == 'Mut' & TET2_6_mut1 == 'Mut'  ~ "Undetermined",
    BCOR_4_insertion1 == 'Undetected' & BCOR_4_insertion2_1 == 'Undetected' &  ASXL1_13_mut == 'Mut' &  RUX1_5_mut == 'Mut' &  TET2_6_mut1 == 'Mut' ~ "Undetermined",
    ASXL1_13_mut == 'Mut' & BCOR_4_insertion1 == 'Mut' & BCOR_4_insertion2_1 %in% c("WT","Undetected")  & RUX1_5_mut == 'Mut' & TET2_6_mut1 == 'Mut'  ~ "TARB1",
    ASXL1_13_mut == 'Mut' & BCOR_4_insertion1 %in% c("WT","Undetected") & BCOR_4_insertion2_1 == 'Mut' & RUX1_5_mut == 'Mut' & TET2_6_mut1 == 'Mut'  ~ "TARB2",
    ASXL1_13_mut == 'Mut' & BCOR_4_insertion1 %in% c("WT","Undetected") & BCOR_4_insertion2_1 %in% c("WT","Undetected") & RUX1_5_mut == 'Mut' & TET2_6_mut1 == 'Mut'  ~ "TAR",
    ASXL1_13_mut == 'Mut' & BCOR_4_insertion1 == 'Mut' & BCOR_4_insertion2_1 %in% c("WT","Undetected") & RUX1_5_mut %in% c("WT","Undetected") & TET2_6_mut1 == 'Mut'  ~ "TARB1",
    ASXL1_13_mut == 'Mut' & BCOR_4_insertion1 %in% c("WT","Undetected") & BCOR_4_insertion2_1 == 'Mut' & RUX1_5_mut %in% c("WT","Undetected") & TET2_6_mut1 == 'Mut'  ~ "TARB2",
    ASXL1_13_mut %in% c("WT","Undetected") & BCOR_4_insertion1 == 'Mut' & BCOR_4_insertion2_1 %in% c("WT","Undetected") & RUX1_5_mut == 'Mut' & TET2_6_mut1 == 'Mut'  ~ "TARB1",
    ASXL1_13_mut %in% c("WT","Undetected") & BCOR_4_insertion1 %in% c("WT","Undetected") & BCOR_4_insertion2_1 == 'Mut' & RUX1_5_mut == 'Mut' & TET2_6_mut1 == 'Mut'  ~ "TARB2",
    ASXL1_13_mut %in% c("WT","Undetected") & BCOR_4_insertion1 == 'Mut' & BCOR_4_insertion2_1 %in% c("WT","Undetected") & RUX1_5_mut %in% c("WT","Undetected") & TET2_6_mut1 == 'Mut'  ~ "TARB1",
    ASXL1_13_mut %in% c("WT","Undetected") & BCOR_4_insertion1 %in% c("WT","Undetected") & BCOR_4_insertion2_1 == 'Mut' & RUX1_5_mut %in% c("WT","Undetected") & TET2_6_mut1 == 'Mut'  ~ "TARB2",
    ASXL1_13_mut == 'Mut' & BCOR_4_insertion1 %in% c("WT","Undetected") & BCOR_4_insertion2_1 %in% c("WT","Undetected") & RUX1_5_mut %in% c("WT","Undetected") & TET2_6_mut1 == 'Mut'  ~ "TA",
    ASXL1_13_mut %in% c("WT","Undetected") & BCOR_4_insertion1 %in% c("WT","Undetected") & BCOR_4_insertion2_1 %in% c("WT","Undetected") & RUX1_5_mut %in% c("WT","Undetected") & TET2_6_mut1 == 'Mut'  ~ "T",
    ASXL1_13_mut == 'Mut' & BCOR_4_insertion1 == 'Mut' & BCOR_4_insertion2_1 %in% c("WT","Undetected") & RUX1_5_mut == 'Mut' & TET2_6_mut1 %in% c("WT","Undetected")  ~ "TARB1",
    ASXL1_13_mut == 'Mut' & BCOR_4_insertion1 %in% c("WT","Undetected") & BCOR_4_insertion2_1 == 'Mut' & RUX1_5_mut == 'Mut' & TET2_6_mut1 %in% c("WT","Undetected")  ~ "TARB2",
    ASXL1_13_mut == 'Mut' & BCOR_4_insertion1 %in% c("WT","Undetected") & BCOR_4_insertion2_1 %in% c("WT","Undetected") & RUX1_5_mut == 'Mut' & TET2_6_mut1 %in% c("WT","Undetected")  ~ "TAR",
    ASXL1_13_mut == 'Mut' & BCOR_4_insertion1 == 'Mut' & BCOR_4_insertion2_1 %in% c("WT","Undetected") & RUX1_5_mut %in% c("WT","Undetected") & TET2_6_mut1 %in% c("WT","Undetected")  ~ "TARB1",
    ASXL1_13_mut == 'Mut' & BCOR_4_insertion1 %in% c("WT","Undetected") & BCOR_4_insertion2_1 == 'Mut' & RUX1_5_mut %in% c("WT","Undetected") & TET2_6_mut1 %in% c("WT","Undetected")  ~ "TARB2",
    ASXL1_13_mut == 'WT' & BCOR_4_insertion1 == 'WT' & BCOR_4_insertion2_1 == 'WT' & RUX1_5_mut == 'WT' & TET2_6_mut1 == 'WT'  ~ "no-somatic-mut",
    ASXL1_13_mut %in% c("WT","Undetected") & BCOR_4_insertion1 %in% c("WT","Undetected") & BCOR_4_insertion2_1 %in% c("WT","Undetected") & RUX1_5_mut == 'Mut' & TET2_6_mut1 == 'Mut'  ~ "TR",
    ASXL1_13_mut == 'Mut' & BCOR_4_insertion1 == 'WT' & BCOR_4_insertion2_1 == 'WT' & RUX1_5_mut == 'WT' & TET2_6_mut1 %in% c("WT","Undetected")  ~ "TA",
    ASXL1_13_mut == 'WT' & BCOR_4_insertion1 == 'WT' & BCOR_4_insertion2_1 == 'WT' & RUX1_5_mut == 'Mut' & TET2_6_mut1 == 'WT'  ~ "TAR",
    ASXL1_13_mut == 'WT' & BCOR_4_insertion1 == 'WT' & BCOR_4_insertion2_1 == 'Undetected' & RUX1_5_mut == 'WT' & TET2_6_mut1 == 'WT'  ~ "no-somatic-mut",
    ASXL1_13_mut == 'WT' & BCOR_4_insertion1 == 'Undetected' & BCOR_4_insertion2_1 == 'WT' & RUX1_5_mut == 'WT' & TET2_6_mut1 == 'WT'  ~ "no-somatic-mut",
    TRUE ~ "Undetermined"))

write_csv(genotypes, 'FINAL_genotypes.csv') # This contains clonal information that will be integrated with scATAC data

# For numbers of cells in a clone
clone_count <- genotypes %>% 
  group_by(clone) %>%
  summarize(count=n())
clone_count %>% write_csv('BRI2349_clone_count_reanalysis.csv')
