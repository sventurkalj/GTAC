#Differential analyses between clones

library(ArchR)
library(parallel)
library(readr)
library(dplyr)
library(chromVAR)
library(chromVARmotifs)
library(ggrepel)
library(Seurat)
library(EnhancedVolcano)

addArchRThreads(threads = 2) 

addArchRGenome("hg38")
peaks_on_clones <- loadArchRProject(path = "./", force = FALSE, showLogo = TRUE)

#Differential gene activity between clones within a cluster
#Subdivide a cluster into multiple groups depending on clone

peaks_on_clones$which_LMPP <- case_when(peaks_on_clones$cell_type == 'LMPP-like' & peaks_on_clones$clone == 'TAR' ~ "LMPP-runx1",
                                         peaks_on_clones$cell_type == 'LMPP-like' & peaks_on_clones$clone == 'TARB1' ~ "LMPP-bcor1",
                                         peaks_on_clones$cell_type == 'LMPP-like' & peaks_on_clones$clone == 'TARB2' ~ "LMPP-bcor2",
                                         TRUE ~ 'other'
)

gene_active_lmpp_bcor2 <- getMarkerFeatures(
  ArchRProj = peaks_on_clones, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "which_LMPP",
  useGroups = 'LMPP-runx1',
  bgdGroups = 'LMPP-bcor2',
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(gene_active_lmpp_bcor2, cutOff = "FDR <= 0.05")
m <- markerList@listData[["LMPP1-runx1"]]@listData %>% 
  as.data.frame() #here you have a table for only significant genes

markerList <- getMarkers(gene_active_lmpp_bcor2, cutOff = "FDR <= 0.5")

name <- markerList@listData[["LMPP1-runx1"]]@listData[["name"]]
value <- markerList@listData[["LMPP1-runx1"]]@listData[["Log2FC"]]
fdr <- markerList@listData[["LMPP1-runx1"]]@listData[["FDR"]]
more_open <- data.frame(name,value,fdr)
more_open$minus_log_10 <- -log10(fdr)

more_open <- more_open%>%
  mutate(updown = case_when(
    value >= 0.25 & fdr <= 0.05 ~ 'UP',
    value < 0.25 & value > -0.25 ~ 'unchanged',
    fdr > 0.05 ~ 'unchanged',
    TRUE ~ 'DOWN'
  ))

p <- EnhancedVolcano(more_open,
                     lab = name,
                     x = 'value',
                     y = 'fdr',
                     selectLab = markerMotifs,
                     xlab = 'Genes TAR vs TARB2 LMPP',
                     pCutoff = 0.05,
                     FCcutoff = 0.25,
                     pointSize = 1.5,
                     labSize = 1,
                     col=c('grey', 'grey', 'grey', 'red3'),
                     colAlpha = 1,
                     drawConnectors = TRUE,
                     boxedLabels = TRUE,
                     maxoverlapsConnectors = 100,
                     xlim = c(-2.1, 2.1),
                     ylim = c(0,18)
) 
plotPDF(p, name = 'genes LMPP.pdf', width = 8, height = 8, addDOC = FALSE)

#Visualize gene openness on a track within the same cluster, between two clones

marker <- "INPP4B"

p <- plotBrowserTrack(
  ArchRProj = peaks_on_clones, 
  groupBy = "which_LMPP1",
  geneSymbol = marker, 
  upstream = 600000,
  downstream = 200000,
  ylim = c(0,1)
)

#GSEA on differentially accessible genes between clones

markerList <- getMarkers(gene_active_lmpp_bcor2, cutOff = "FDR <= 1")
DE_genes <- markerList@listData[["LMPP1-runx1"]]@listData %>% 
  as.data.frame() 

DE_genes$ranking <- (-log10(DE_genes$FDR))*sign(DE_genes$Log2FC)

# Order list by the ranking
DE_genes<-DE_genes[order(DE_genes$ranking,decreasing = T),]

prerank.genes<- DE_genes$ranking
names(prerank.genes)<-DE_genes$name

# Get list of gene sets from gmt file

gene_sets<- fgsea::gmtPathways("c2.cgp.v2022.1.Hs.symbols.gmt") 
gene_sets<- fgsea::gmtPathways("c2.cp.v2022.1.Hs.symbols.gmt")
gene_sets<- fgsea::gmtPathways("c3.tft.v2022.1.Hs.symbols.gmt")
gene_sets<- fgsea::gmtPathways("h.all.v2022.1.Hs.symbols.gmt")
gene_sets<- fgsea::gmtPathways("c5.go.v2022.1.Hs.symbols.gmt")

# Run analysis with fgsea

set.seed(42)
fgseaRes <- fgsea::fgseaMultilevel(pathways = gene_sets, 
                                   stats = prerank.genes,
                                   minSize=5,
                                   maxSize=1000,
                                   eps = 1e-50,
                                   nPermSimple=10000) %>%
  filter(padj < 0.05)

p <- fgsea::plotEnrichment(gene_sets[["EPPERT_CE_HSC_LSC"]],
                           stats= prerank.genes)
plotPDF(p, name="EPPERT_CE_HSC_LSC_TARBvsTAR_LMPP.pdf")

#Pairwise test for differentially accessible peaks

markerTest <- getMarkerFeatures(
  ArchRProj = peaks_on_clones, 
  useMatrix = "PeakMatrix",
  groupBy = "which_LMPP",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "LMPP-runx1",
  bgdGroups = "LMPP-bcor2"
)

markerList <- getMarkers(markerTest, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5")
m <- markerList@listData[["LMPP1-runx1"]]@listData %>% 
  as.data.frame()
write_csv(m, 'Revision_peaks_0.25LMPP1_TARvsTARB2.csv')

#For pairwise differential openness visualization of motif accessibility, using chromVAR values

diffMotif <- getMarkerFeatures(
  ArchRProj = peaks_on_clones, 
  testMethod = "wilcoxon",
  useGroups = "LMPP-runx1",
  bgdGroups = "LMPP-bcor2",
  binarize = FALSE,
  useMatrix = "MotifMatrix",
  groupBy = "which_LMPP",
  useSeqnames="z"
)

markerList <- getMarkers(diffMotif, cutOff = "FDR <= 0.5")

name <- markerList@listData[["LMPP1-runx1"]]@listData[["name"]]
value <- markerList@listData[["LMPP1-runx1"]]@listData[["MeanDiff"]]
fdr <- markerList@listData[["LMPP1-runx1"]]@listData[["FDR"]]
more_open <- data.frame(name,value,fdr)
more_open$minus_log_10 <- -log10(fdr)

more_open <- more_open%>%
  mutate(updown = case_when(
    value >= 0.25 & fdr <= 0.05 ~ 'UP',
    value < 0.25 & value > -0.25 ~ 'unchanged',
    fdr > 0.05 ~ 'unchanged',
    TRUE ~ 'DOWN'
  ))

motifs_TARB2 <-c('DNMT1', 'SMAD5', 'RELA', 'REL', 'NFKB1', 'NFKB2', 'VENTX', 'NANOG', 'FOS', 'JUN', 'FOSB', 'JUNB', 'JUND', 'IRF4', 'CEBPA')
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))

p <- EnhancedVolcano(more_open,
                     lab = name,
                     x = 'value',
                     y = 'fdr',
                     selectLab = markerMotifs,
                     xlab = 'Genes TAR vs TARB2 LMPP',
                     pCutoff = 0.05,
                     FCcutoff = 0.25,
                     pointSize = 1.5,
                     labSize = 1,
                     col=c('grey', 'grey', 'grey', 'red3'),
                     colAlpha = 1,
                     drawConnectors = TRUE,
                     boxedLabels = TRUE,
                     maxoverlapsConnectors = 100,
                     xlim = c(-2.1, 2.1),
                     ylim = c(0,18)
) 
plotPDF(p, name = 'motifs LMPP.pdf', width = 8, height = 8, addDOC = FALSE)

#TF footprinting - clone-specific TF footprints

motifPositions <- getPositions(peaks_on_clones)

motifs <- c("ETV6", "ERG", "WT1", "PBX3", "MEF2D", "FLI1", "SMAD5", "FOXP2", "HOXA4", "HOXD1", "FOXO1", "FOXO3", "MEF2C", "TCF7", "NFE2", "REL", "NFKB1", "NFKB2", "RELA", "HOXD8", "HOXB3", "HOXB5", "NANOG", "VENTX", "RELB", "HOXB2", "HOXB7", "HOXB8", "IRF4", "BATF", "ID1", "ID2", "BCL11A", "BCL11B", "MYCL1", "MYC", "CEBPA", "CEBPD")
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))

peaks_on_clones <- addGroupCoverages(ArchRProj = peaks_on_clones, groupBy = "which_LMPP", force=TRUE)

seFoot <- getFootprints(
  ArchRProj = peaks_on_clones, 
  positions = motifPositions[markerMotifs], 
  groupBy = "which_LMPP",
  useGroups = c("LMPP-runx1", "LMPP-bcor2")
)

plotFootprints(
  seFoot = seFoot,
  ArchRProj = peaks_on_clones, 
  normMethod = "Subtract",
  plotName = "Subtracted_LMPP_footprints.pdf",
  addDOC = FALSE,
  smoothWindow = 5
)

#Trajectory analysis 1: clone-specific trajectories
#These cell labellings will be used to generate clone-specific trajectories for three lineages within leukemia 
#Note that: LMPP1 = LMPP; LMPP2 = LM progenitor-like; GMP1: pDC/NKP-like; GMP2: Mono/Neut progenitor-like; PreB = B progenitor like

peaks_on_clones$TARB2_trajectory_cells <- case_when(peaks_on_clones$which_LMPP1 == 'LMPP1-bcor2'  ~ "LMPP1-bcor2",
                                                    peaks_on_clones$which_LMPP2 == 'LMPP2-bcor2'  ~ "LMPP2-bcor2", 
                                                    peaks_on_clones$which_GMP1 == 'GMP1-bcor2'  ~ "GMP1-bcor2",
                                                    peaks_on_clones$which_GMP2 == 'GMP2-bcor2'  ~ "GMP2-bcor2",
                                                    peaks_on_clones$which_preB == 'PreB-bcor2'  ~ "PreB-bcor2",
                                                    TRUE ~ 'other'
)

peaks_on_clones$TARB1_trajectory_cells <- case_when(peaks_on_clones$which_LMPP1 == 'LMPP1-bcor1'  ~ "LMPP1-bcor1",
                                                    peaks_on_clones$which_LMPP2 == 'LMPP2-bcor1'  ~ "LMPP2-bcor1", 
                                                    peaks_on_clones$which_GMP1 == 'GMP1-bcor1'  ~ "GMP1-bcor1",
                                                    peaks_on_clones$which_GMP2 == 'GMP2-bcor1'  ~ "GMP2-bcor1",
                                                    peaks_on_clones$which_preB == 'PreB-bcor1'  ~ "PreB-bcor1",
                                                    TRUE ~ 'other'
)

peaks_on_clones$TAR_trajectory_cells <- case_when(peaks_on_clones$which_LMPP1 == 'LMPP1-runx1'  ~ "LMPP1-runx1",
                                                  peaks_on_clones$which_LMPP2 == 'LMPP2-runx1'  ~ "LMPP2-runx1", 
                                                  peaks_on_clones$which_GMP1 == 'GMP1-runx1'  ~ "GMP1-runx1",
                                                  peaks_on_clones$which_GMP2 == 'GMP2-runx1'  ~ "GMP2-runx1",
                                                  peaks_on_clones$which_preB == 'PreB-runx1'  ~ "PreB-runx1",
                                                  TRUE ~ 'other'
)

#We will show an example of clone-specific analysis within the pDC-NKP-like lineage

#Get TAR-only trajectory towards pDC/NKP-like

trajectory <- c("LMPP1-runx1", "LMPP2-runx1", "GMP1-runx1")

peaks_on_clones <- addTrajectory(
  ArchRProj = peaks_on_clones, 
  name = "pDC_NKP_TAR", 
  groupBy = "TAR_trajectory_cells",
  trajectory = trajectory, 
  embedding = "UMAP", 
  force = TRUE
)

head(peaks_on_clones$pDC_NKP_TAR[!is.na(peaks_on_clones$pDC_NKP_TAR)])
p <- plotTrajectory(peaks_on_clones, trajectory = "pDC_NKP_TAR", colorBy = "cellColData", name = "pDC_NKP_TAR")
p[[1]]
plotPDF(p, name = "TAR_pDC_NKP_trajectory.pdf", addDOC = FALSE, width = 5, height = 5)

#Get TARB2-only trajectory towards pDC/NKP-like

trajectory <- c("LMPP1-bcor2", "LMPP2-bcor2", "GMP1-bcor2")
trajectory

peaks_on_clones <- addTrajectory(
  ArchRProj = peaks_on_clones, 
  name = "pDC_NKP_TARB2", 
  groupBy = "TARB2_trajectory_cells",
  trajectory = trajectory, 
  embedding = "UMAP", 
  force = TRUE
)

head(peaks_on_clones$pDC_NKP_TARB2[!is.na(peaks_on_clones$pDC_NKP_TARB2)])
p <- plotTrajectory(peaks_on_clones, trajectory = "pDC_NKP_TARB2", colorBy = "cellColData", name = "pDC_NKP_TARB2")
p[[1]]
plotPDF(p, name = "TARB2_pDC_NKP_trajectory.pdf", addDOC = FALSE, width = 5, height = 5)

#Get TARB1-only trajectory towards pDC/NKP-like

trajectory <- c("LMPP1-bcor1", "LMPP2-bcor1", "GMP1-bcor1")
trajectory

peaks_on_clones <- addTrajectory(
  ArchRProj = peaks_on_clones, 
  name = "pDC_NKP_TARB1", 
  groupBy = "TARB1_trajectory_cells",
  trajectory = trajectory, 
  embedding = "UMAP", 
  force = TRUE
)

head(peaks_on_clones$pDC_NKP_TARB1[!is.na(peaks_on_clones$pDC_NKP_TARB1)])
p <- plotTrajectory(peaks_on_clones, trajectory = "pDC_NKP_TARB1", colorBy = "cellColData", name = "pDC_NKP_TARB1")
p[[1]]
plotPDF(p, name = "TARB1_pDC_NKP_trajectory.pdf", addDOC = FALSE, width = 5, height = 5)

#Infer gene  dynamics in trajectory

p1 <- plotTrajectory(peaks_on_clones, trajectory = "pDC_NKP_TAR", colorBy = "GeneScoreMatrix", name = "IRF8", continuousSet = "horizonExtra", plotAs = "points", size = 2.3)
p1[[1]]
plotPDF(p1, name = 'IRF8_pDC_TAR.pdf', addDOC = FALSE, width = 5, height = 5) #Repeat this for other two clone-specific trajectories

#Infer motif dynamics in trajectory

p1 <- plotTrajectory(peaks_on_clones, trajectory = "pDC_NKP_TAR", colorBy = "MotifMatrix", name = "z:BCL11A_194", continuousSet = "horizonExtra", plotAs = "points", size = 2.3)
p1[[1]]
plotPDF(p1, name = 'dc_tar_bcl11a.pdf', addDOC = FALSE, width = 5, height = 5)

#Infer global motif dynamics in mutant trajectories in an unbiased way

trajMM  <- getTrajectory(ArchRProj = peaks_on_clones, name = "pDC_NKP_TAR", useMatrix = "MotifMatrix", log2Norm = FALSE)
p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"), labelTop = 200)
plotPDF(p1, name = "pDC_NKP_TAR_motifs.pdf", addDOC = FALSE, width = 6, height = 14)

#Infer global gene dynamics

trajGSM <- getTrajectory(ArchRProj = peaks_on_clones, name = "pDC_NKP_TAR", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
p2 <- trajectoryHeatmap(trajGSM,  pal = paletteContinuous(set = "horizonExtra"), labelTop = 100)
plotPDF(p2, name = "pDC_NKP_TAR.pdf", addDOC = FALSE, width = 6, height = 14)

#Split clone-specific trajectories into five quantiles to study how genomic elements vary through differentiation in clone-specific manner

peaks_on_clones$pDC_NKP_TAR_TARB2_quantiles <- case_when(peaks_on_clones$pDC_NKP_TAR > 0 & peaks_on_clones$pDC_NKP_TAR < 20 ~ "WT1",
                                              peaks_on_clones$pDC_NKP_TAR > 20 & peaks_on_clones$pDC_NKP_TAR < 40 ~ "WT2", 
                                              peaks_on_clones$pDC_NKP_TAR > 40 & peaks_on_clones$pDC_NKP_TAR < 60 ~ "WT3", 
                                              peaks_on_clones$pDC_NKP_TAR > 60 & peaks_on_clones$pDC_NKP_TAR < 80 ~ "WT4",
                                              peaks_on_clones$pDC_NKP_TAR > 80 & peaks_on_clones$pDC_NKP_TAR <= 100 ~ "WT5",
                                              peaks_on_clones$pDC_NKP_TARB2 > 0 & peaks_on_clones$pDC_NKP_TARB2 < 20 ~ "B1",
                                              peaks_on_clones$pDC_NKP_TARB2 > 20 & peaks_on_clones$pDC_NKP_TARB2 < 40 ~ "B2", 
                                              peaks_on_clones$pDC_NKP_TARB2 > 40 & peaks_on_clones$pDC_NKP_TARB2 < 60 ~ "B3", 
                                              peaks_on_clones$pDC_NKP_TARB2 > 60 & peaks_on_clones$pDC_NKP_TARB2 < 80 ~ "B4",
                                              peaks_on_clones$pDC_NKP_TARB2 > 80 & peaks_on_clones$pDC_NKP_TARB2 <= 100 ~ "B5",
                                              TRUE ~ 'other'
)

#Now you can run ArchR Browser and generate tracks using this column

#Trajectory analysis 2: generate non clone-specific trajectories - example for pDC/NKP lineage

trajectoryNK <- c("LMPP-like-1", "LMPP-like-2", "GMP-like-1")

peaks_on_clones <- addTrajectory(
  ArchRProj = peaks_on_clones, 
  name = "pDC_NKP_overall_trajectory", 
  groupBy = "cell_type",
  trajectory = trajectoryNK, 
  embedding = "UMAP", 
  force = TRUE
)

head(peaks_on_clones$pDC_NKP_overall_trajectory[!is.na(peaks_on_clones$pDC_NKP_overall_trajectory)])
p <- plotTrajectory(peaks_on_clones, trajectory = "pDC_NKP_overall_trajectory", colorBy = "cellColData", name = "pDC_NKP_overall_trajectory", continuousSet = "blueYellow")
plotPDF(p, name = "pDC_traj.pdf", addDOC = FALSE, width = 5, height = 5)

#Plot clonal presence across pseudotime

archR_meta <- as.data.frame(peaks_on_clones@cellColData)
clones <- archR_meta[,20] #Take clonal information
pseudotime <- archR_meta[,54] #pDC/NKP pseudotime
pseudotime_pDC <- data.frame(clones,pseudotime) %>%
  write_csv("total_nk_traj.csv") 

library(cowplot)
theme_set(theme_cowplot())

p <- ggplot(pseudotime_pDC, aes(x = pseudotime))+
  geom_density(aes(colour = clones, fill = clones),alpha = 0.3)+
  theme_cowplot() + labs(x = 'Clones_pDC_Pseudotime')
plotPDF(p, name = "Clones_pDC_Pseudotime.pdf", addDOC = FALSE, width = 6.5, height = 6.5)

# Understading dynamics of genes and motifs through full trajectories

# Get motif matrix from project
matrix <- getMatrixFromProject(
  ArchRProj = peaks_on_clones,
  useMatrix = "MotifMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = FALSE,
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
)

names <- matrix@colData@rownames %>% #Get cell names
  as.data.frame 

m <- matrix@assays
m <- m@data@listData[["z"]]
m <- as.matrix(m) %>%
  t() %>%
  as.data.frame() #This is to get the matrix of z deviations x cells

archR_meta <- tibble::rownames_to_column(archR_meta, "cell_ID")
m <- tibble::rownames_to_column(m, "cell_ID")
all_motifs <- left_join(archR_meta,m, by='cell_ID')

#This plots clone-specific motif accessibility through pseudotime
p <- ggplot(all_motifs, aes(x = pDC_NKP_overall_trajectory, y = NFKB1_719))+
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), aes(color = clone))+
  theme_cowplot() + labs(x = 'pDC pseudotime', y = 'NFKB1_719 motif deviation')
plotPDF(p, name = "NFKB1_719 pDC gam.pdf", addDOC = FALSE, width = 6.5, height = 6.5)

# Get gene matrix from project
matrix <- getMatrixFromProject(
  ArchRProj = peaks_on_clones,
  useMatrix = "GeneScoreMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = FALSE,
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
)

# Some work to extract per-cell gene scores for each gene
m <- matrix@assays@data@listData[["GeneScoreMatrix"]] %>%
  as.matrix() %>%
  as.data.frame()
names <- matrix@elementMetadata@listData[["name"]]  %>%
  as.data.frame() 

names <- tibble::rownames_to_column(names, "ID")
m <- tibble::rownames_to_column(m, "ID")

genes_cells <- left_join(names,m, by='ID')
genes_cells$ID <- NULL
names(genes_cells)[names(genes_cells) == '.'] <- 'gene'

rownames(genes_cells) <- genes_cells[,1]
genes_cells$gene <- NULL

genes_cells1 <- t(genes_cells) %>%
  as.matrix() %>%
  as.data.frame()
genes_cells1 <- tibble::rownames_to_column(genes_cells1, "ID")

archR_meta <- as.data.frame(peaks_on_clones@cellColData)
archR_meta <- tibble::rownames_to_column(archR_meta, "ID")
total_gene_cells <- left_join(archR_meta, genes_cells1, by = 'ID') # Generates a table containing metadata and gene scores for all genes for each cell

# Use the table to plot clone-specific gene activity score for any gene across the given trajectory

p <- ggplot(total_genes_cells, aes(x = pDC_NKP_overall_trajectory, y = INPP4B))+
  geom_smooth(aes(colour = clone)) +
  theme_cowplot() + labs(x = 'pDC_NKP_pseudotime', y = 'INPP4B gene activity')
plotPDF(p, name = "INPP4B_pDC_NKP_pseudotime.pdf", addDOC = FALSE, width = 6.5, height = 6.5)

# Correlate accessibility of TF motifs across clusters

library(ggpubr)

motif <- all_motifs[c("cell_type","GFI1_220","SMAD5_866")] %>% # Take any two motifs of interest and cell type information
  subset(cell_type == 'LMPP-like-1') # Select cell type of interest

p <- ggplot(motif, aes(x=GFI1_220, y=SMAD5_866)) + 
  geom_point()+
  geom_smooth(method=lm) +
  theme_cowplot() +
  stat_cor(method = "spearman", label.x = 1, label.y = 5)
ggsave("GFI1_220_smad5_lmpp_corr.pdf", device = "pdf", width = 7)

# Compare motif accessibility for genotypes

motif <- all_motifs[c("which_LMPP1","SMAD5_866")] %>% # Take any two motifs of interest and cell type information
  subset(which_LMPP1 == 'LMPP1-runx1' | which_LMPP1 == 'LMPP1-bcor2') # Select cell type of interest

mycol <- c("green", "orange")
plot1 <- ggplot(motif, aes(x=which_LMPP1, y=SMAD5_866)) +
  geom_boxplot(fill = mycol) +
  xlab("LMPP genotype")+
  ylab("SMAD5 motif accessibility")+
  stat_compare_means()+
  theme_cowplot()
ggsave("smad5_lmpp.pdf", device = 'pdf')
