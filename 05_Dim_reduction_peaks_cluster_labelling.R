#ArchR analysis for dimensionality reduction, clustering, peak calling and cluster assignment

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
AML <- loadArchRProject(path = "./", force = FALSE, showLogo = TRUE)

#Dimensionality reduction - run LSI

AML <- addIterativeLSI(
  ArchRProj = AML,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2), 
    sampleCells = 4000, 
    n.start = 10
  ), 
  varFeatures = 30000, 
  dimsToUse = 1:30,
  force = TRUE
)

#We had one sample so did not run Harmony. However, if you have multiple samples, it is worth running it now to use it for integration

AML <- addHarmony(
  ArchRProj = AML,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample"
)

#Clustering using Seurat's methods - this is deterministic, so every time the same input will generate the same output. One can play with resolution here, depending on how granular cell populations need to be.

AML <- addClusters(
  input = AML,
  reducedDims = "IterativeLSI", #Alternatively this can be Harmony for example
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8,
  force = TRUE
)

#Generate UMAP

AML <- addUMAP(
  ArchRProj = AML, 
  reducedDims = "IterativeLSI", #Alternatively this can be Harmony for example
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force = TRUE
)

#Label UMAP by cell clusters
p1 <- plotEmbedding(ArchRProj = AML,
                    colorBy = "cellColData", 
                    name = "Clusters", 
                    embedding = "UMAP", 
                    size = 1, 
                    labelMeans = FALSE
)
plotPDF(p1, name = "AML_clusters.pdf", addDOC = FALSE, width = 5, height = 5)

#Label UMAP by clone, to see how genotypes are distributed on the embedding
p2 <- plotEmbedding(ArchRProj = AML,
                    colorBy = "cellColData", 
                    name = "clone", 
                    embedding = "UMAP", 
                    size = 1, 
                    labelMeans = FALSE
)
plotPDF(p1, name = "AML_clones.pdf", addDOC = FALSE, width = 5, height = 5)

#Label UMAP by plate, to confirm the absence of plate-related batch effects
p3 <- plotEmbedding(ArchRProj = AML,
                    colorBy = "cellColData", 
                    name = "plate", 
                    embedding = "UMAP", 
                    size = 1, 
                    labelMeans = FALSE
)
plotPDF(p3, name = "AML_plates.pdf", addDOC = FALSE, width = 5, height = 5)

#At this point perform a Fisher's exact test to see whether there is significant association of some clones with some clusters

#Get marker genes for clusters or clones based on GENE ACTIVITY SCORE 

markersGS <- getMarkerFeatures(
  ArchRProj = AML, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Cluster",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 0.5") #You can experiment with the Log2FC thresholds here depending on how well the populations are resolved

markers <- as.data.frame(markerList) %>%
  write_csv('marker_genes_0.5log2fc.csv') #Export the list

#Generate heatmap with marker genes for each cluster

markerGenes  <- c(
  "CD34", 'CD38', 'MEIS1', 'HLF', 'AVP','MECOM', 'HOXB8', 'GATA2', 'HOXA9', 'ELANE', 'CTSG', 'CEBPA',  #Early Progenitor
  "GATA1", 'LST1', 'TAL1', 'CLIP3', #Erythroid
  "PAX5", "MS4A1", "EBF1", "MME", 'ETS1', 'TNFRSF4','TNFRSF9', 'CLEC6A', 'PRDM16', 'TPO', 'IRF8', 'TLR4',    #B-Cell Trajectory
  "CD14", "CEBPB", 'PBX1', 'DNTT','BLNK', 'TPSB2', 'PTPN14', 'SPIB', 'SIGLEC5', 'PLD1', 'SIGLEC10', 'IRF4', 'CEBPA', 'CEBPB',  #Monocytes
  "IRF8", 'MPO', 'HBB', 'NFIB', 'PF4', 'MLLT3', 'MPL', 'VWF', 'GYPA', 'GATA1', 'CD3D', 'EOMES', 'LEF1', 'CD8A', 'CD3E', 
  "CD3D", "CD8A", "TBX21", "IL7R", 'CD4', 'BCL11B', 'BLNK', 'CD79B', 'CD79A', 'MS4A1', 'IRF4', 'FLT3', 'AIF1', 'GP9', 'GZMB', 'PRF1', 'CD9', 'PLEK', 'ITGA2B', 'CD33', 'RNASE3'  #TCells 
)

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.05 & Log2FC >= 0.5", 
  labelMarkers = markerGenes,
  transpose = TRUE, clusterCols = F 
)
p <- ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(p, name='0.5_cluster_Heatmap_gene_activity.pdf', width = 10)

#Visualize markers on embedding - MAGIC used for signal smoothing

markerGenes  <- c(
  "MN1", "CD19", 'ERG', "BAALC", "CD22", 'CD79A'
) #Any genes you are interested to visualize

AML <- addImputeWeights(AML) #These will be used to smooth the signal, but it is also good to check with the non-smoothed graphs

p <- plotEmbedding(
  ArchRProj = AML, 
  colorBy = "GeneScoreMatrix", 
  name = marker, 
  embedding = "UMAP",
  plotAs = 'points',
  quantCut = c(0.01, 0.95),
  size = 2.3
)
ggsave('markers_without_smoothing.pdf', device = "pdf", height = 5, width = 5)

p <- plotEmbedding(
  ArchRProj = AML, 
  colorBy = "GeneScoreMatrix", 
  name = marker, 
  embedding = "UMAP",
  plotAs = 'points',
  quantCut = c(0.01, 0.95),
  imputeWeights = getImputeWeights(peaks_on_clones),
  size = 2.3
)
ggsave('markers_with_smoothing.pdf', device = "pdf", height = 5, width = 5)

#Track plotting with ArchR browser

marker <- "CEBPD"
p <- plotBrowserTrack(
  ArchRProj = AML, 
  groupBy = "cluster",
  geneSymbol = marker, 
  upstream = 600000,
  downstream = 200000,
  ylim = c(0,1) #To visualize a specific gene

)

coordinates <- GRanges(seqnames='chrX', ranges = IRanges(start=2691310, width=42363))
p <- plotBrowserTrack(
  ArchRProj = AML, 
  groupBy = "cluster",
  region = coordinates, 
  ylim = c(0,1)
) #To visualize a specific genomic area
plotPDF(p, "CD99_track.pdf")

#Visualize tracks by clone
marker <- "CEBPD"
p <- plotBrowserTrack(
  ArchRProj = AML, 
  groupBy = "clone",
  geneSymbol = marker, 
  upstream = 600000,
  downstream = 200000,
  ylim = c(0,1) #To visualize a specific gene in a clone-specific manner
  
)

#To launch the browser
ArchRBrowser(AML)

#GSEA on differentially accessible genes - to understand which hematopoietic signatures are enriched in which clusters

AML$is_cluster_1 <- case_when(AML$cluster == '1'  ~ "yes",
                                            TRUE ~ 'no' #You can do this for every cluster to compare them against everything else
)

gene_active <- getMarkerFeatures(
  ArchRProj = AML, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "is_cluster_1",
  useGroups = 'yes',
  bgdGroups = 'no',
  maxCells = 2000,
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(gene_active, cutOff = "FDR <= 1")
DE_genes <- markerList@listData[["yes"]]@listData %>% 
  as.data.frame() 

DE_genes$ranking <- (-log10(DE_genes$FDR))*sign(DE_genes$Log2FC)

# Order list by the ranking
DE_genes<-DE_genes[order(DE_genes$ranking,decreasing = T),]

prerank.genes<- DE_genes$ranking
names(prerank.genes)<-DE_genes$name

# Get list of gene sets from gmt file - here you can input any set of signatures

gene_sets<- fgsea::gmtPathways("human.hematopoiesis.signature.genes.curated.v6.gmt") 

# Run analysis with fgsea

set.seed(42)
fgseaRes <- fgsea::fgseaMultilevel(pathways = gene_sets, 
                                   stats = prerank.genes,
                                   minSize=5,
                                   maxSize=1000,
                                   eps = 1e-50,
                                   nPermSimple=10000) %>%
  filter(padj < 0.05)

bubble <- bubble %>%
  mutate(value = -log10(padj)) %>%
  arrange(desc(NES)) %>%
  mutate(pathway = fct_reorder(pathway, NES))
p1 <- ggplot(bubble, aes(x=NES, y = pathway, size = value, color = NES )) +
  geom_point() + 
  theme_cowplot()+
  colorspace::scale_color_continuous_divergingx(palette = 'RdBu', mid = 0, limits = c(-3.5,3.5), rev = T)
ggsave("enrichment cluster 1.pdf",p1, device = "pdf", width = 9)

#Make pseudobulk replicates - this is key to call peaks. We made our pseudobulk replicates based on clones, but one can do them based on clusters.
#If you have more clones than clusters, it might be worth making them by clone. However, usually you have more clusters.
#The best thing is to do both, call peaks using both and see what result is more sensible. 

AML <- addGroupCoverages(ArchRProj = AML, groupBy = "clone")

pathToMacs2 <- findMacs2() #transfer the project and code to the server and use macs2 there
pathToMacs2 <- '/path_to_macs2'

AML <- addReproduciblePeakSet(
  ArchRProj = AML, 
  groupBy = "clone", 
  pathToMacs2 = pathToMacs2
)

getPeakSet(AML)
AML <- addPeakMatrix(AML) 

peaks_on_clones <- AML

#Get marker peaks for clusters. Later these will be used to identify marker TF motifs for each cluster

markersPeaks_clusters <- getMarkerFeatures(
  ArchRProj = peaks_on_clones, 
  useMatrix = "PeakMatrix", 
  groupBy = "cluster",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

#Get marker peaks for clones. Later these will be used to identify marker TF motifs for each clone.

markersPeaks_clones <- getMarkerFeatures(
  ArchRProj = peaks_on_clones, 
  useMatrix = "PeakMatrix", 
  groupBy = "clone",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

#Motif analysis in differentially accessible peaks

peaks_on_clones <- addMotifAnnotations(ArchRProj = peaks_on_clones, motifSet = "cisbp", name = "Motif")

#Set markerTest on lines 277 and 287; 

motifs_clusters <- peakAnnoEnrichment(   
  seMarker = markersPeaks_clusters,
  ArchRProj = peaks_on_clones,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.05 & Log2FC >= 0.5" #can experiment with the Log2FC value
)

motif_clones <- peakAnnoEnrichment(
  seMarker = markersPeaks_clones,
  ArchRProj = peaks_on_clones,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.05 & Log2FC >= 0.5"
)

heatmapEM <- plotEnrichHeatmap(motifs_clusters, n = 30, transpose = TRUE, rastr = FALSE)
m <- ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(m, name='marker_motifs.pdf', width = 12)

#Custom enrichments - overlay marker peaks with published hematopoietic bulk ATAC-seq data

peaks_on_clones <- addArchRAnnotations(ArchRProj = peaks_on_clones, collection = "ATAC")

enrichATAC <- peakAnnoEnrichment(
  seMarker = markersPeaks_clusters,
  ArchRProj = peaks_on_clones,
  peakAnnotation = "ATAC",
  cutOff = "FDR <= 0.05 & Log2FC >= 0.25"
)

heatmapATAC <- plotEnrichHeatmap(enrichEncode, n = 15, transpose = TRUE, rastr = FALSE)
m <- ComplexHeatmap::draw(heatmapEncode, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(m, name='0.25_clusters_ATAC anno.pdf')

#ChromVAR to get motif deviations

if("Motif" %ni% names(peaks_on_clones@peakAnnotation)){
  peaks_on_clones <- addMotifAnnotations(ArchRProj = peaks_on_clones, motifSet = "cisbp", name = "Motif")
}

peaks_on_clones <- addBgdPeaks(peaks_on_clones) #background peaks for deviation computing

peaks_on_clones <- addDeviationsMatrix(
  ArchRProj = peaks_on_clones, 
  peakAnnotation = "Motif",
  force = TRUE
)

#For UMAP visualization of deviations

motifs <- c("EOMES", "HOXA9", 'CEBPA','LEF1', 'GATA1', 'SPI1', 'EBF1', 'PAX5', 'ETV6', 'RUNX1', 'FOS', 'TBX21', 'MEF2C', 'MEF2D') #Set here TFs likely to distinguish different populations
markerMotifs <- getFeatures(peaks_on_clones, select = paste(motif, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs <- grep("z:", markerMotifs, value = TRUE)

p <- plotEmbedding(
  ArchRProj = peaks_on_clones, 
  colorBy = "MotifMatrix", 
  name = sort(markerMotifs), 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(peaks_on_clones),
  size = 2.3,
  plotAs = 'points'
)
plotPDF(p, name = 'motifs.pdf')

p2 <- lapply(p, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
plot <- do.call(cowplot::plot_grid, c(list(ncol = 3),p2))
plotPDF(plot, name='chromVAR_UMAP.pdf')

#Assigning cluster names: we have enough evidence now to call clusters

peaks_on_clones$cell_type <- case_when(peaks_on_clones$clusters == '1'  ~ "LMPP-like",
                                            peaks_on_clones$clusters == '2'  ~ "LM progenitor-like", 
                                            peaks_on_clones$clusters == '3'  ~ "DC/NK progenitor-like",
                                            peaks_on_clones$clusters == '4'  ~ "Mono/Neut progenitor-like",
                                            peaks_on_clones$clusters == '5'  ~ "B progenitor-like",
                                            peaks_on_clones$clusters == '6'  ~ "HSC/MPP",
                                            peaks_on_clones$clusters == '7'  ~ "Mk/EP",
                                            peaks_on_clones$clusters == '8'  ~ "T/NK cells",
                                            TRUE ~ 'other'
)

#Module scores
peaks_on_clones <- addImputeWeights(peaks_on_clones) #If not done already

LSC <- list(leukaemia = c('SETBP1', 'HLF', 'GIMAP7', 'ABCC2', 'HOPX', 'TMEM200A', 'PCDHGC3', 'FCMR', 'RBPMS', 'EFCC1', 'LTB', 'EBF3', 'GUCY1A1', 'SH3BP5', 'SLC37A3', 'GIMAP6', 'CD34', 'VNN1', 'MMRN1', 'SLC38A1', 'GIMAP2', 'BIRC3', 'LGALSL', 'EVI2A', 'GSAP', 'FBXO21', 'MEF2C', 'ICAM1', 'HECA'),
            dummy = c('HLF', 'HOXA9', 'MECOM', 'AVP')
)

peaks_on_clones <- addModuleScore(peaks_on_clones, features = LSC, useMatrix = "GeneScoreMatrix")
p1 <- plotEmbedding(peaks_on_clones, name="leukaemia", imputeWeights = getImputeWeights(peaks_on_clones))
plotPDF(ggAlignPlots(p1,p2,draw=F,type="h"))