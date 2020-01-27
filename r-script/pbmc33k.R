#Note : This is not a full tutorial but an annotated command list
#For a more detailed tutorial and intro to Seurat, please see our tutorial of the 3k dataset

library(Seurat)
library(Matrix)
pbmc33k.data <- Read10X("~/Projects/datasets/pbmc33k/filtered_gene_bc_matrices/hg19/")
pbmc33k  <- new("seurat", raw.data = pbmc33k.data)

#setup setting do.scale and do.center to F - this means that we will NOT scale genes by default (to speed things up)
pbmc33k <- Setup(pbmc33k, min.cells = 3, min.genes = 200, project = "10X_PBMC", do.scale = F, do.center = F, names.field = 2,
	names.delim = "\\-")
mito.genes <- grep("^MT-", rownames(pbmc33k@data), value = T)
percent.mito <- colSums(expm1(pbmc33k@data[mito.genes, ])) / colSums(expm1(pbmc33k@data))

#AddMetaData adds columns to object@data.info, and is a great place to stash QC stats
pbmc33k <- AddMetaData(pbmc33k, percent.mito, "percent.mito")
VlnPlot(pbmc33k, c("nGene", "nUMI", "percent.mito"), nCol = 3)


#We filter out cells that have unique gene counts over 2,500 and under 500, and > 5% mitochondrial percentage
pbmc33k <- SubsetData(pbmc33k, subset.name = "nGene", accept.high = 2500)
pbmc33k <- SubsetData(pbmc33k, subset.name = "percent.mito", accept.high = 0.05)
pbmc33k <- SubsetData(pbmc33k, subset.name = "nGene", accept.low = 500)

#choose gene outliers on mean-variability plot
pbmc33k <- MeanVarPlot(pbmc33k, x.low.cutoff = 0, y.cutoff = 0.8)
length(pbmc33k@var.genes)

#Perform negative-binomial regression on the variable genes, this sets their value in pbmc33k@scale.data, which is used for PCA/clustering
#We only do this on the variable genes to save time, but you can do this genome-wide
#We treat mitochondrial percentage, batch, and nUMI as confounding variables, 

#you can save the object at any time to save results, and can restore it back in using load()
save(pbmc33k, file = "~/Projects/datasets/pbmc33k/pbmc33k_tutorial.Robj")

pbmc33k <- RegressOut(pbmc33k,latent.vars = c("percent.mito", "orig.ident", "nUMI"), genes.regress = pbmc33k@var.genes, model.use = "negbinom")

#Run PCA with the IRLBA package (iteratively computes the top dimensions, dramatic increase in speed since we are throwing away most PCs anyway)
pbmc33k <- PCAFast(pbmc33k, pc.genes = pbmc33k@var.genes, pcs.compute = 40, pcs.print = 30)

PCElbowPlot(pbmc33k, num.pc = 40)
PrintPCA(pbmc33k,pcs.print = 1:36)
PCHeatmap(pbmc33k,pc.use = 1:12,100)
PCHeatmap(pbmc33k,pc.use = 13:24,100)
PCHeatmap(pbmc33k,pc.use = 25:36,100)

#select 25 PCs for downstream analysis
pbmc33k <- RunTSNE(pbmc33k, dims.use = 1:25, do.fast = T)

#save.SNN means that you can easily re-run with different resolution values. Here we run with a few different res v
pbmc33k <- FindClusters(pbmc33k ,pc.use = 1:25, resolution = seq(2,4,0.5), save.SNN = T, do.sparse = T)

#pbmc33k_20=FindClusters(pbmc33k,pc.use = 1:25,resolution = seq(0.6,4,0.1),save.SNN = T,do.sparse = T,k.param = 20)
save(pbmc33k,file = "~/Projects/datasets/pbmc33k/pbmc33k_nbinom_final.Robj")

#Explore results with a resolution value of 2
#We see nice clusters, but the data is slightly under-clustered (for example CD1C+ and CD141+ DCs are merged together)
pbmc33k <- SetAllIdent(pbmc33k, id = "res.2")
TSNEPlot(pbmc33k,do.label = T)


#In the section below, we explore an overclustering combined with post-hoc merging strategy that can help discover weaker splits in the data
#This section applies cutoffs in an admittedly supervised way, and is something we are actively working to improve
#As a more conservative appraoch, ignore the section below, and use a resolution value of 2 as shown above.

#We can bump the resolution up to call more clusters, but this slightly over-clusters the data
pbmc33k <- SetAllIdent(pbmc33k, id = "res.4")
TSNEPlot(pbmc33k,do.label = T)

#In fact there is no 'perfect' value of resolution, we always either under or over-cluster the data
#This is because of dramatically different cluster sizes, and is known as the 'multi-resolution' problem in graph-based clustering
#One solution is to slightly over-cluster the data, and then perform a post-hoc merging step, where transcriptionally indistinguishable clusters are merged back together
#As a test for merging, we use the Out-of-bag error (OOBE) from a random forest classifier, but you could also set a cutoff for # of differentially expressed genes

pbmc33k <- SetAllIdent(pbmc33k, id = "res.4")

#Build a classification hierarchy that places transcriptionally similar clusters adjacent on a tree
pbmc33k <- BuildClusterTree(pbmc33k, do.reorder = T, reorder.numeric = T)

#calculate the classification error for left/right cells on each branch of the tree
#sort internal nodes based on OOBE. For nodes with high OOBE, we cannot accurately tell the left/right children apart based on random forests, so the clusters may need to be merged
node.scores <- AssessNodes(pbmc33k)
node.scores[order(node.scores$oobe,decreasing = T),] -> node.scores

#choose the first eight splits to merge )
nodes.merge=node.scores[1:8,]
nodes.to.merge=sort(nodes.merge$node)
pbmc33k.merged <- pbmc33k
for (n in nodes.to.merge){
  pbmc33k.merged <- MergeNode(pbmc33k.merged, n)
}

#rebuild the classification hierarchy using all genes for interpretation (you can also try with variable genes as the default)
pbmc33k.merged <- BuildClusterTree(pbmc33k.merged, do.reorder = T, reorder.numeric = T,genes.use = rownames(pbmc33k.merged@data))
TSNEPlot(pbmc33k.merged, do.label = T)

#Rebuild the tree to update the leaf identities (upcoming version will update this automatically)
pbmc33 <- BuildClusterTree(pbmc33k.merged, genes.use = rownames(pbmc33k.merged@data))

#Other commands that might be useful. See function documentation for more complete descriptions
PlotClusterTree(pbmc33k.merged)

#color TSNE based on a hierarchical split
ColorTSNESplit(pbmc33k.merged, node = 33)
ColorTSNESplit(pbmc33k.merged, node = 31,color1 = "red",color3="blue")

#Visualize canonical markers
FeaturePlot(o4, c("MS4A1", "GNLY","CD3E","CD8A","LYZ","PF4"),cols.use = c("lightgrey","blue"),nCol = 3)

#Find myeloid/lymphoid markers, speed up by requiring a 25% difference in detection rate between cell populations
#Can set test.use="negbinom" to apply our new negative binomial DE test
FindMarkersNode(pbmc33.merged, node = 31, min.diff.pct = 0.25)

#Dot plot visualization
DotPlot(pbmc33k.merged,c("CD3E","CD8B","SELL","GNLY","GZMA","GZMB","GZMH","GZMK","PRF1","NKG7","XCL2","FGFBP2","KLRC1"),cols.use = myPalette(low="lightgrey",high = "blue"),cex.use = 2)

#For easy exploration, create a smaller object with a max of 100 cells per cluster
#Some functions in beta to display heatmap markers of each split in the tree
pbmc33k.downsampled <- SubsetData(pbmc33k.merged, max.cells.per.ident = 100)
markers <- FindAllMarkersNode(pbmc33k.downsampled, node = 31, max.cells.per.ident = 100)
HeatmapNode(pbmc33k.downsampled, node = 34, marker.list = markers, max.genes = 10)
FindMarkers(pbmc33k.downsampled, ident.1 = 17,  only.pos = T)


#rename cluster IDs
new_ids=c("Megakaryocyte","Mono_CD16+","Mono_CD16+_C1qa","Mono_Apobec3b","Mono_Apobec3a","Mono_CD14+_Antiviral","Mono_CD14+","Mono_CD14+_Inflam","DC_CD141+",
          "MonoT_Doublet","DC_CD1C+","NK_Prss57","NK","NK_XCL1","Multiplet","CD8Effector_GZMH","CD8Effector_GZMK","BCell_A","BCell_B","BCell_C","CD34+",
          "T_DoubleNeg","CD8_Memory","CD8_Naive","CD4_Naive","CD4_Memory_Antiviral","CD4_Memory_TIGIT","CD4_Memory","pDC","Plasma")

current.cluster.ids=1:30
pbmc33k.merged@ident <- plyr::mapvalues(pbmc33k.merged@ident, from = current.cluster.ids, to = new_ids)


#Find markers for CD1C vs CD141 DCs, using the negative binomial test
#Can speed up if desired  by setting max.cells.per.ident or min.diff.pct
dc.markers=FindMarkers(pbmc33k.merged,"DC_CD1C+","DC_CD141+",test.use = "negbinom")



