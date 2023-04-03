#Input scATAC data into ArchR project and merge with genotyping data
#### Load the relevant libraries, create arrow files and ArchR project  ####
library(ArchR)
library(parallel)
library(readr)
library(dplyr)
library(chromVAR)
library(chromVARmotifs)
library(ggrepel)
library(Seurat)

addArchRThreads(threads = 2) 

addArchRGenome("hg38")

# Format fragment file so that it is compatible with ArchR
inputFiles <- "fragment_file.tsv.gz"
reformatFragmentFiles(
  fragmentFiles = inputFiles,
  checkChrPrefix = getArchRChrPrefix()
)

#Create Arrow Files and set initial QC thershold - note that we used a rather high minFrags thershold
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = "AML",
  maxFrags = 120000,
  minTSS = 8, 
  minFrags = 10000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  excludeChr = c("chrM", "chrY"),
  force=TRUE
)

#Infer doublets with ArchR's method 
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, 
  knnMethod = "UMAP", 
  LSIMethod = 1
)

#Create ArchR project
AML <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "AML_GTAC_analysis",
  copyArrows = FALSE #keeps an original copy of arrow files for later usage.
)
archR_meta <- as.data.frame(AML@cellColData) %>%
  write_csv("AML_initial_meta.csv") #Read ArchR metadata and observe the initial QC results

#### Adding genotype data into the ArchR project ####
geno_metadata <- read.csv('FINAL_genotypes.csv') #This is the file with clonal information for each cell generated in the previous section
archR_meta$ID <- rownames(archR_meta) #Convert rownames to column to use it for table joining

meta1 <- archR_meta %>% 
  separate(col = ID, sep = '#', into = c('a','ID')) #This is to get rid of the '#' symbol in cell names provided by ArchR

#Convert the cell ID in geno_metadata to exactly match the one provided by ArchR - this will obviously depend on naming strtagy for genotyping data
geno_metadata <- geno_metadata %>% 
  separate(col = X, sep = '_', into = c('a','b','c','d','e','f')) %>%
  unite('ID', b:f, sep='_') 

total_metadata <- full_join(archR_meta, sven_metadata, by = "ID") #Join the two tables via cell IDs but KEEPING the exact order present in ArchR metadata
total_metadata <- total_metadata[-c(3538:3820),] #Substract the final rows which are not present in ArchR metadata - these are cells excluded by initial QC
rownames(total_metadata) <- total_metadata$ID

#Now you have clonal information for each cell in the ArchR project. We want to input this data into the ArchR project so that it stays saved for future use. 

clone <- total_metadata$clone   #This will generate a vector from total_metadata column, RESPECTING the exact cell order from ArchR metadata
AML$clone <- clone

#Clonal info is now in the ArchR project. 
#Equally to clone info, we use this way to input any metadata of interest into ArchR. You can repeat this to add info about plate, sample etc. 
#Discard cells from project for which clone is unknown. For analysis, we only used cells for which the clone was known 

idxSample <- BiocGenerics::which(AML$clone %in% NA)
to_discard <- AML$cellNames[-idxSample]
AML <- AML[to_discard, ]

#We also discarded a few cells for which clone was TR
idxSample <- BiocGenerics::which(AML$clone %in% "TR")
to_discard <- AML$cellNames[-idxSample]
AML <- AML1[to_discard, ]

#Finally we also discarded cells coming from the WT bone marrow control
idxSample <- BiocGenerics::which(AML$which_sample %in% "NOC153") 
to_discard <- AML$cellNames[-idxSample]
AML <- AML[to_discard, ]

#Plotting QC metrics
tss_fragments <- getCellColData(AML, select = c("log10(nFrags)", "TSSEnrichment"))

p <- ggPoint(
  x = tss_fragments[,1], 
  y = tss_fragments[,2], 
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), log10(500000)),
  ylim = c(0, 26)
) + geom_hline(yintercept = 8, lty = "dashed") + geom_vline(xintercept = 4, lty = "dashed") 
plotPDF(p, name = "AML-TSS-vs-Frags.pdf", addDOC = FALSE)

#Plot metrics by plate to get an idea whether there could be a bath effect or difference in quality

p1 <- plotGroups(
  ArchRProj = AML, 
  groupBy = "plate", 
  colorBy = "cellColData", 
  name = "DoubletEnrichment",
  maxCells = 10000,
  plotAs = "violin"
)

p2 <- plotGroups(
  ArchRProj = AML, 
  groupBy = "plate", 
  colorBy = "cellColData", 
  name = "nFrags",
  maxCells = 10000,
  plotAs = "violin"
)

p3 <- plotGroups(
  ArchRProj = AML, 
  groupBy = "plate", 
  colorBy = "cellColData", 
  name = "NucleosomeRatio",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)

#Plotting fragment size distribution and TSS

p4 <- plotFragmentSizes(ArchRProj = AML, threads = 4)
plotPDF(p4, name='fragment_size.pdf')

p5 <- plotTSSEnrichment(ArchRProj = AML, threads=4)
plotPDF(p5, name='TSS.pdf')

#Filter out doublets (did with a modified function written by Jeff Granja on github error thread, use filter 1)

AML <- filterDoublets(AML)
saveArchRProject(ArchRProj = AML, outputDirectory = "./", load = FALSE)