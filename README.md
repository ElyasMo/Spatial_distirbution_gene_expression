# Spatial_distirbution_gene_expression
Here I have developed a pipeline to monitor the spatial distribution of gene expression in various Seurat objects based on four different thresholds. 

**The R markdown output**
---
title: "TMEM119-spatial distribution of expression"
author: "Elyas Mohammadi"
date: "April 2, 2021"
output:
  html_document:
---

```{r}

library(Seurat)
library(sctransform)

#####Reading in the 10xVisium objects


setwd("D:/Poland/PHD/spatial/Second_set/HE_sample/")
Directory="D:/Poland/PHD/spatial/Second_set/HE_sample/"
HE_sample_2<-Load10X_Spatial(Directory,filename = "filtered_feature_bc_matrix.h5",
                             assay = "Spatial",
                             slice = "HEsmpl2",
                             filter.matrix = TRUE,
                             to.upper = FALSE)

setwd("D:/Poland/PHD/spatial/Second_set/HE_control/")
Directory="D:/Poland/PHD/spatial/Second_set/HE_control/"
HE_control_2<-Load10X_Spatial(Directory,filename = "filtered_feature_bc_matrix.h5",
                              assay = "Spatial",
                              slice = "HEctr2",
                              filter.matrix = TRUE,
                              to.upper = FALSE)

setwd("D:/Poland/PHD/spatial/project/HE/HE_sample_2248/")
Directory="D:/Poland/PHD/spatial/project/HE/HE_sample_2248/"
HE_sample_1<-Load10X_Spatial(Directory,filename = "filtered_feature_bc_matrix.h5",
                             assay = "Spatial",
                             slice = "HEsmpl1",
                             filter.matrix = TRUE,
                             to.upper = FALSE)

setwd("D:/Poland/PHD/spatial/project/HE/HE_control_2339/")
Directory="D:/Poland/PHD/spatial/project/HE/HE_control_2339/"
HE_control_1<-Load10X_Spatial(Directory,filename = "filtered_feature_bc_matrix.h5",
                              assay = "Spatial",
                              slice = "HEctr1",
                              filter.matrix = TRUE,
                              to.upper = FALSE)
#####integration
setwd("D:/Poland/PHD/spatial/project/CR/CR_sample_2248/")
Directory="D:/Poland/PHD/spatial/project/CR/CR_sample_2248/"
CR_sample_1<-Load10X_Spatial(Directory,filename = "filtered_feature_bc_matrix.h5",
                             assay = "Spatial",
                             slice = "CRsmpl1",
                             filter.matrix = TRUE,
                             to.upper = FALSE)
setwd("D:/Poland/PHD/spatial/project/CR/CR_sample_2248/")
Directory="D:/Poland/PHD/spatial/project/CR/CR_sample_2248/"
CR_sample_1_density<-Load10X_Spatial(Directory,filename = "filtered_feature_bc_matrix.h5",
                             assay = "Spatial",
                             slice = "CRsmpl1",
                             filter.matrix = TRUE,
                             to.upper = FALSE)
setwd("D:/Poland/PHD/spatial/project/CR/CR_control_2339/")
Directory="D:/Poland/PHD/spatial/project/CR/CR_control_2339/"
CR_control_1<-Load10X_Spatial(Directory,filename = "filtered_feature_bc_matrix.h5",
                              assay = "Spatial",
                              slice = "CRctr1",
                              filter.matrix = TRUE,
                              to.upper = FALSE)

##Background
setwd("D:/Poland/PHD/spatial/Second_set/HE_sample/")
Directory="D:/Poland/PHD/spatial/Second_set/HE_sample/"
HE_sample_2_back<-Load10X_Spatial(Directory,filename = "filtered_feature_bc_matrix.h5",
                             assay = "Spatial",
                             slice = "HEsmpl2",
                             filter.matrix = TRUE,
                             to.upper = FALSE)

setwd("D:/Poland/PHD/spatial/Second_set/HE_control/")
Directory="D:/Poland/PHD/spatial/Second_set/HE_control/"
HE_control_2_back<-Load10X_Spatial(Directory,filename = "filtered_feature_bc_matrix.h5",
                              assay = "Spatial",
                              slice = "HEctr2",
                              filter.matrix = TRUE,
                              to.upper = FALSE)

setwd("D:/Poland/PHD/spatial/project/HE/HE_sample_2248/")
Directory="D:/Poland/PHD/spatial/project/HE/HE_sample_2248/"
HE_sample_1_back<-Load10X_Spatial(Directory,filename = "filtered_feature_bc_matrix.h5",
                             assay = "Spatial",
                             slice = "HEsmpl1",
                             filter.matrix = TRUE,
                             to.upper = FALSE)

setwd("D:/Poland/PHD/spatial/project/HE/HE_control_2339/")
Directory="D:/Poland/PHD/spatial/project/HE/HE_control_2339/"
HE_control_1_back<-Load10X_Spatial(Directory,filename = "filtered_feature_bc_matrix.h5",
                              assay = "Spatial",
                              slice = "HEctr1",
                              filter.matrix = TRUE,
                              to.upper = FALSE)
#####integration
setwd("D:/Poland/PHD/spatial/project/CR/CR_sample_2248/")
Directory="D:/Poland/PHD/spatial/project/CR/CR_sample_2248/"
CR_sample_1_back<-Load10X_Spatial(Directory,filename = "filtered_feature_bc_matrix.h5",
                             assay = "Spatial",
                             slice = "CRsmpl1",
                             filter.matrix = TRUE,
                             to.upper = FALSE)

setwd("D:/Poland/PHD/spatial/project/CR/CR_control_2339/")
Directory="D:/Poland/PHD/spatial/project/CR/CR_control_2339/"
CR_control_1_back<-Load10X_Spatial(Directory,filename = "filtered_feature_bc_matrix.h5",
                              assay = "Spatial",
                              slice = "CRctr1",
                              filter.matrix = TRUE,
                              to.upper = FALSE)


##adding the percentage of MT genes to metadata
CR_control_1[["percent.mt"]] <- PercentageFeatureSet(CR_control_1, pattern = "^MT-")
CR_sample_1[["percent.mt"]] <- PercentageFeatureSet(CR_sample_1, pattern = "^MT-")
CR_sample_1_density[["percent.mt"]] <- PercentageFeatureSet(CR_sample_1_density, pattern = "^MT-")
HE_control_1[["percent.mt"]] <- PercentageFeatureSet(HE_control_1, pattern = "^MT-")
HE_sample_1[["percent.mt"]] <- PercentageFeatureSet(HE_sample_1, pattern = "^MT-")
HE_control_2[["percent.mt"]] <- PercentageFeatureSet(HE_control_2, pattern = "^MT-")
HE_sample_2[["percent.mt"]] <- PercentageFeatureSet(HE_sample_2, pattern = "^MT-")

CR_control_1_back[["percent.mt"]] <- PercentageFeatureSet(CR_control_1_back, pattern = "^MT-")
CR_sample_1_back[["percent.mt"]] <- PercentageFeatureSet(CR_sample_1_back, pattern = "^MT-")
HE_control_1_back[["percent.mt"]] <- PercentageFeatureSet(HE_control_1_back, pattern = "^MT-")
HE_sample_1_back[["percent.mt"]] <- PercentageFeatureSet(HE_sample_1_back, pattern = "^MT-")
HE_control_2_back[["percent.mt"]] <- PercentageFeatureSet(HE_control_2_back, pattern = "^MT-")
HE_sample_2_back[["percent.mt"]] <- PercentageFeatureSet(HE_sample_2_back, pattern = "^MT-")



##Filtering the metadata basd on various criteria
CR_control_1 <- subset(CR_control_1, subset = nFeature_Spatial > 200 & nFeature_Spatial < 7000 & percent.mt < 15)
CR_sample_1 <- subset(CR_sample_1, subset = nFeature_Spatial > 200 & nFeature_Spatial < 7000 & percent.mt < 15 )
CR_sample_1_density <- subset(CR_sample_1_density, subset = nFeature_Spatial > 200 & nFeature_Spatial < 7000 & percent.mt < 15 )
HE_control_1 <- subset(HE_control_1, subset = nFeature_Spatial > 200 & nFeature_Spatial < 7000 & percent.mt < 15)
HE_sample_1 <- subset(HE_sample_1, subset = nFeature_Spatial > 200 & nFeature_Spatial < 7000 & percent.mt < 15)
HE_control_2 <- subset(HE_control_2, subset = nFeature_Spatial > 200 & nFeature_Spatial < 7000 & percent.mt < 15)
HE_sample_2 <- subset(HE_sample_2, subset = nFeature_Spatial > 200 & nFeature_Spatial < 7000 & percent.mt < 15)


CR_control_1_back <- subset(CR_control_1_back, subset = nFeature_Spatial > 200 & nFeature_Spatial < 7000 & percent.mt < 15)
CR_sample_1_back <- subset(CR_sample_1_back, subset = nFeature_Spatial > 200 & nFeature_Spatial < 7000 & percent.mt < 15 )
HE_control_1_back <- subset(HE_control_1_back, subset = nFeature_Spatial > 200 & nFeature_Spatial < 7000 & percent.mt < 15)
HE_sample_1_back <- subset(HE_sample_1_back, subset = nFeature_Spatial > 200 & nFeature_Spatial < 7000 & percent.mt < 15)
HE_control_2_back <- subset(HE_control_2_back, subset = nFeature_Spatial > 200 & nFeature_Spatial < 7000 & percent.mt < 15)
HE_sample_2_back <- subset(HE_sample_2_back, subset = nFeature_Spatial > 200 & nFeature_Spatial < 7000 & percent.mt < 15)




#SCtransform function-Normalization
HE_sample_2<-SCTransform(HE_sample_2, assay = "Spatial", verbose = FALSE)
HE_sample_2<- FindVariableFeatures(HE_sample_2, selection.method = "vst", 
                                   nfeatures = 3000, verbose = FALSE)

HE_control_2<-SCTransform(HE_control_2, assay = "Spatial", verbose = FALSE)
HE_control_2<- FindVariableFeatures(HE_control_2, selection.method = "vst", 
                                    nfeatures = 3000, verbose = FALSE)

CR_sample_1 <- SCTransform(CR_sample_1, assay = "Spatial", verbose = FALSE)
CR_sample_1<- FindVariableFeatures(CR_sample_1, selection.method = "vst", 
                                   nfeatures = 3000, verbose = FALSE)

CR_sample_1_density <- SCTransform(CR_sample_1_density, assay = "Spatial", verbose = FALSE)
CR_sample_1_density<- FindVariableFeatures(CR_sample_1_density, selection.method = "vst", 
                                   nfeatures = 3000, verbose = FALSE)

CR_control_1<-SCTransform(CR_control_1, assay = "Spatial", verbose = FALSE)
CR_control_1<- FindVariableFeatures(CR_control_1, selection.method = "vst", 
                                    nfeatures = 3000, verbose = FALSE)

HE_sample_1<-SCTransform(HE_sample_1, assay = "Spatial", verbose = FALSE)
HE_sample_1<- FindVariableFeatures(HE_sample_1, selection.method = "vst", 
                                   nfeatures = 3000, verbose = FALSE)
HE_control_1<-SCTransform(HE_control_1, assay = "Spatial", verbose = FALSE)
HE_control_1<- FindVariableFeatures(HE_control_1, selection.method = "vst", 
                                    nfeatures = 3000, verbose = FALSE)



HE_sample_2_back<-SCTransform(HE_sample_2_back, assay = "Spatial", verbose = FALSE)
HE_sample_2_back<- FindVariableFeatures(HE_sample_2_back, selection.method = "vst", 
                                   nfeatures = 3000, verbose = FALSE)

HE_control_2_back<-SCTransform(HE_control_2_back, assay = "Spatial", verbose = FALSE)
HE_control_2_back<- FindVariableFeatures(HE_control_2_back, selection.method = "vst", 
                                    nfeatures = 3000, verbose = FALSE)

CR_sample_1_back <- SCTransform(CR_sample_1_back, assay = "Spatial", verbose = FALSE)
CR_sample_1_back <- FindVariableFeatures(CR_sample_1_back, selection.method = "vst", 
                                   nfeatures = 3000, verbose = FALSE)


CR_control_1_back<-SCTransform(CR_control_1, assay = "Spatial", verbose = FALSE)
CR_control_1_back<- FindVariableFeatures(CR_control_1, selection.method = "vst", 
                                    nfeatures = 3000, verbose = FALSE)

HE_sample_1_back<-SCTransform(HE_sample_1, assay = "Spatial", verbose = FALSE)
HE_sample_1_back<- FindVariableFeatures(HE_sample_1, selection.method = "vst", 
                                   nfeatures = 3000, verbose = FALSE)
HE_control_1_back<-SCTransform(HE_control_1_back, assay = "Spatial", verbose = FALSE)
HE_control_1_back<- FindVariableFeatures(HE_control_1_back, selection.method = "vst", 
                                    nfeatures = 3000, verbose = FALSE)

```


```{r, warning = FALSE, message = FALSE, cache=TRUE, fig.show='hold', out.width="100%", fig.fullwidth = TRUE}
#Insertin the desired gene or list of genes in addition to four diferent threasholds.
Genes <- ("TMEM119")


##Calculating the normalized percentage of desired genes in each spot
AD_gene <- match(Genes,rownames(CR_sample_1@assays$SCT))
AD_genes <- rownames(CR_sample_1@assays$SCT)[AD_gene]
AD_genes <- AD_genes[!is.na(AD_genes)]
CR_sample_1[["percent.AD"]]<-PercentageFeatureSet(CR_sample_1,features=AD_genes)

AD_gene <- match(Genes,rownames(CR_control_1@assays$SCT))
AD_genes <- rownames(CR_control_1@assays$SCT)[AD_gene]
AD_genes <- AD_genes[!is.na(AD_genes)]
CR_control_1[["percent.AD"]]<-PercentageFeatureSet(CR_control_1,features=AD_genes)

AD_gene <- match(Genes,rownames(HE_sample_1@assays$SCT))
AD_genes <- rownames(HE_sample_1@assays$SCT)[AD_gene]
AD_genes <- AD_genes[!is.na(AD_genes)]
HE_sample_1[["percent.AD"]]<-PercentageFeatureSet(HE_sample_1,features=AD_genes)

AD_gene <- match(Genes,rownames(HE_control_1@assays$SCT))
AD_genes <- rownames(HE_control_1@assays$SCT)[AD_gene]
AD_genes <- AD_genes[!is.na(AD_genes)]
HE_control_1[["percent.AD"]]<-PercentageFeatureSet(HE_control_1,features=AD_genes)

AD_gene <- match(Genes,rownames(HE_sample_2@assays$SCT))
AD_genes <- rownames(HE_sample_2@assays$SCT)[AD_gene]
AD_genes <- AD_genes[!is.na(AD_genes)]
HE_sample_2[["percent.AD"]]<-PercentageFeatureSet(HE_sample_2,features=AD_genes)

AD_gene <- match(Genes,rownames(HE_control_2@assays$SCT))
AD_genes <- rownames(HE_control_2@assays$SCT)[AD_gene]
AD_genes <- AD_genes[!is.na(AD_genes)]
HE_control_2[["percent.AD"]]<-PercentageFeatureSet(HE_control_2,features=AD_genes)


##A good estimation for the threasholds
##Box plot may help to decide about the threasholds

boxplot(CR_sample_1@meta.data$percent.AD, CR_control_1@meta.data$percent.AD, HE_sample_1@meta.data$percent.AD, HE_control_1@meta.data$percent.AD, HE_sample_2@meta.data$percent.AD, HE_control_2@meta.data$percent.AD , main = "Multiple boxplots for comparision of % of desired gene", names=c("CR_SP1", "CR_CT1", "HE_SP1", "HE_CT1", "HE_SP2", "HE_CT2"), col = c( "dark blue","orange","red", "yellow", "light blue", "purple"), las = 2)

```
![Boxplot1- Helps to determine threasholds in the next chunck](https://github.com/ElyasMo/Spatial_distirbution_gene_expression/blob/main/Boxplot1.png)
**The threasholds can be added in this section**
```{r, warning = FALSE, message = FALSE, cache=TRUE, fig.show='hold', out.width="100%", fig.fullwidth = TRUE}
THR1 = 0
THR2 = 0.025
THR3 = 0.04
THR4 = 0.05
```

**Adding the cordinates**
The metadata is filtered and we would like to see the locations which have the accumulated spots with high percentage of determined genes

```{r , warning = FALSE, message = FALSE, cache=TRUE}
##inserting the cordinates from spatial directory of spaceranger output
setwd("D:/Poland/PHD/spatial/Second_set/cordinations/")
CR_CT1_spatial <- read.csv("D:/Poland/PHD/spatial/Second_set/cordinations/CR_CT1_spatial.csv", row.names=1)
CR_SP1_spatial <- read.csv("D:/Poland/PHD/spatial/Second_set/cordinations/CR_SP1_spatial..csv", row.names=1)
CR_SP1_density <- read.csv("D:/Poland/PHD/spatial/Second_set/cordinations/CR_SP1_spatial_density.csv", row.names=1)
HE_CT1_spatial <- read.csv("D:/Poland/PHD/spatial/Second_set/cordinations/HE_CT1_spatial.csv", row.names=1)
HE_CT2_spatial <- read.csv("D:/Poland/PHD/spatial/Second_set/cordinations/HE_CT2_spatial.csv", row.names=1)
HE_SP1_spatial <- read.csv("D:/Poland/PHD/spatial/Second_set/cordinations/HE_SP1_spatial.csv", row.names=1)
HE_SP2_spatial <- read.csv("D:/Poland/PHD/spatial/Second_set/cordinations/HE_SP2_spatial.csv", row.names=1)

##Adding cordinates to the metadata
CR_CT1_meta <- cbind(CR_control_1@meta.data[intersect(rownames(CR_control_1@meta.data), rownames(CR_CT1_spatial)),], CR_CT1_spatial[intersect(rownames(CR_control_1@meta.data), rownames(CR_CT1_spatial)),])

CR_SP1_meta <-cbind(CR_sample_1@meta.data[intersect(rownames(CR_sample_1@meta.data), rownames(CR_SP1_spatial)),], CR_SP1_spatial[intersect(rownames(CR_sample_1@meta.data), rownames(CR_SP1_spatial)),])

CR_SP1_density_meta <-cbind(CR_sample_1_density@meta.data[intersect(rownames(CR_sample_1_density@meta.data), rownames(CR_SP1_density)),], CR_SP1_density[intersect(rownames(CR_sample_1_density@meta.data), rownames(CR_SP1_density)),])

CR_SP1_meta <-cbind(CR_sample_1@meta.data[intersect(rownames(CR_sample_1@meta.data), rownames(CR_SP1_spatial)),], CR_SP1_spatial[intersect(rownames(CR_sample_1@meta.data), rownames(CR_SP1_spatial)),])

HE_CT1_meta <-cbind(HE_control_1@meta.data[intersect(rownames(HE_control_1@meta.data), rownames(HE_CT1_spatial)),], HE_CT1_spatial[intersect(rownames(HE_control_1@meta.data), rownames(HE_CT1_spatial)),])

HE_CT2_meta <- cbind(HE_control_2@meta.data[intersect(rownames(HE_control_2@meta.data), rownames(HE_CT2_spatial)),], HE_CT2_spatial[intersect(rownames(HE_control_2@meta.data), rownames(HE_CT2_spatial)),])

HE_SP1_meta <- cbind(HE_sample_1@meta.data[intersect(rownames(HE_sample_1@meta.data), rownames(HE_SP1_spatial)),], HE_SP1_spatial[intersect(rownames(HE_sample_1@meta.data), rownames(HE_SP1_spatial)),])

HE_SP2_meta <- cbind(HE_sample_2@meta.data[intersect(rownames(HE_sample_2@meta.data), rownames(HE_SP2_spatial)),],  HE_SP2_spatial[intersect(rownames(HE_sample_2@meta.data), rownames(HE_SP2_spatial)),])


##Background
CR_CT1_meta_back <- cbind(CR_control_1_back@meta.data[intersect(rownames(CR_control_1_back@meta.data), rownames(CR_CT1_spatial)),], CR_CT1_spatial[intersect(rownames(CR_control_1_back@meta.data), rownames(CR_CT1_spatial)),])

CR_SP1_meta_back <-cbind(CR_sample_1_back@meta.data[intersect(rownames(CR_sample_1_back@meta.data), rownames(CR_SP1_spatial)),], CR_SP1_spatial[intersect(rownames(CR_sample_1_back@meta.data), rownames(CR_SP1_spatial)),])

CR_SP1_meta_back <-cbind(CR_sample_1_back@meta.data[intersect(rownames(CR_sample_1_back@meta.data), rownames(CR_SP1_spatial)),], CR_SP1_spatial[intersect(rownames(CR_sample_1_back@meta.data), rownames(CR_SP1_spatial)),])

HE_CT1_meta_back <-cbind(HE_control_1_back@meta.data[intersect(rownames(HE_control_1_back@meta.data), rownames(HE_CT1_spatial)),], HE_CT1_spatial[intersect(rownames(HE_control_1_back@meta.data), rownames(HE_CT1_spatial)),])

HE_CT2_meta_back <- cbind(HE_control_2_back@meta.data[intersect(rownames(HE_control_2_back@meta.data), rownames(HE_CT2_spatial)),], HE_CT2_spatial[intersect(rownames(HE_control_2_back@meta.data), rownames(HE_CT2_spatial)),])

HE_SP1_meta_back <- cbind(HE_sample_1_back@meta.data[intersect(rownames(HE_sample_1_back@meta.data), rownames(HE_SP1_spatial)),], HE_SP1_spatial[intersect(rownames(HE_sample_1_back@meta.data), rownames(HE_SP1_spatial)),])

HE_SP2_meta_back <- cbind(HE_sample_2_back@meta.data[intersect(rownames(HE_sample_2_back@meta.data), rownames(HE_SP2_spatial)),],  HE_SP2_spatial[intersect(rownames(HE_sample_2_back@meta.data), rownames(HE_SP2_spatial)),])
```

**Modifying the metadata for manual grouping in order to illustrate the amyloid spots**
```{r, warning = FALSE, message = FALSE, cache=TRUE}
##Clustering the CR sample 1 to manually define the amyloid spots
CR_sample_1 <- RunPCA(CR_sample_1, verbose = FALSE)
CR_sample_1 <- RunUMAP(CR_sample_1, dims = 1:30, verbose = FALSE)
CR_sample_1 <- FindNeighbors(CR_sample_1, dims = 1:30, verbose = FALSE)
CR_sample_1 <- FindClusters(CR_sample_1, verbose = FALSE, resolution = 0.3)

clusters <- read.csv("D:/Poland/PHD/spatial/Second_set/cordinations/merged_CRSP1_amy.csv", row.names=1)

CR_sample_1@meta.data['SCT_snn_res.0.3']<-clusters[,6]
CR_sample_1@meta.data['seurat_clusters']<-clusters[,6]

CR_sample_1 <- SetIdent(CR_sample_1, value = CR_sample_1@meta.data$seurat_clusters)


new.cluster.ids <- c('Amyloid','non_Amyloid')
names(new.cluster.ids) <- levels(CR_sample_1)
CR_sample_1 <- RenameIdents(CR_sample_1, new.cluster.ids)


##Here we can see the manual grouping of amyloid spots in first CR sample experiment
SpatialDimPlot(CR_sample_1, label = FALSE, label.size =2, crop = FALSE)
```
![Manually grouped spots](https://github.com/ElyasMo/Spatial_distirbution_gene_expression/blob/main/Manual%20grouping.png)

**Preparing threasholds for plotting**

```{r , warning = FALSE, message = FALSE, cache=TRUE}
library(ggplot2)
library(plyr)
library(viridis)
library(MASS)
library(raster)


##Plotting by kde2d
#kde = kde2d(CR_CT1_meta$Y, -(CR_CT1_meta$X), n = 100)
#contour(kde)

# convert to raster and dataframe
#r = raster(kde)
#df = as.data.frame(r, xy = T)

## group by x value using aggregate() or ddply()
#hd = ddply(df, "x", summarize, max = max(layer))
#agg = aggregate(df$layer, by = list(df$x), FUN = max)

# matching corresponding y values and add to plot
#hd = df[match(hd$max, df$layer),]
#lines(hd$x, hd$y, col = "red", lwd = 2)



#plotting
CR_CT1_meta_THR1 <- subset(CR_CT1_meta, subset =  percent.AD > THR1)
CR_SP1_meta_THR1 <- subset(CR_SP1_meta, subset =  percent.AD > THR1)
HE_CT1_meta_THR1 <- subset(HE_CT1_meta, subset =  percent.AD > THR1)
HE_CT2_meta_THR1 <- subset(HE_CT2_meta, subset =  percent.AD > THR1)
HE_SP1_meta_THR1 <- subset(HE_SP1_meta, subset =  percent.AD > THR1)
HE_SP2_meta_THR1 <- subset(HE_SP2_meta, subset =  percent.AD > THR1)

CR_CT1_meta_THR2 <- subset(CR_CT1_meta, subset =  percent.AD > THR2)
CR_SP1_meta_THR2 <- subset(CR_SP1_meta, subset =  percent.AD > THR2)
HE_CT1_meta_THR2 <- subset(HE_CT1_meta, subset =  percent.AD > THR2)
HE_CT2_meta_THR2 <- subset(HE_CT2_meta, subset =  percent.AD > THR2)
HE_SP1_meta_THR2 <- subset(HE_SP1_meta, subset =  percent.AD > THR2)
HE_SP2_meta_THR2 <- subset(HE_SP2_meta, subset =  percent.AD > THR2)

CR_CT1_meta_THR3 <- subset(CR_CT1_meta, subset =  percent.AD > THR3)
CR_SP1_meta_THR3 <- subset(CR_SP1_meta, subset =  percent.AD > THR3)
HE_CT1_meta_THR3 <- subset(HE_CT1_meta, subset =  percent.AD > THR3)
HE_CT2_meta_THR3 <- subset(HE_CT2_meta, subset =  percent.AD > THR3)
HE_SP1_meta_THR3 <- subset(HE_SP1_meta, subset =  percent.AD > THR3)
HE_SP2_meta_THR3 <- subset(HE_SP2_meta, subset =  percent.AD > THR3)

CR_CT1_meta_THR4 <- subset(CR_CT1_meta, subset =  percent.AD > THR4)
CR_SP1_meta_THR4 <- subset(CR_SP1_meta, subset =  percent.AD > THR4)
HE_CT1_meta_THR4 <- subset(HE_CT1_meta, subset =  percent.AD > THR4)
HE_CT2_meta_THR4 <- subset(HE_CT2_meta, subset =  percent.AD > THR4)
HE_SP1_meta_THR4 <- subset(HE_SP1_meta, subset =  percent.AD > THR4)
HE_SP2_meta_THR4 <- subset(HE_SP2_meta, subset =  percent.AD > THR4)


```

**In this section you can see the plots based on the threashold you determined at the beginning (THR threasholds)**

```{r , warning = FALSE, message = FALSE, cache=TRUE, fig.show='hold', out.width="25%", fig.fullwidth = TRUE}

a <- c("CR_CT1 > THR1, nr spots:","CR_CT1 > THR2, nr spots:","CR_CT1 > THR3, nr spots:","CR_CT1 > THR4, nr spots")
b=1
for (x in c("CR_CT1_meta_THR1","CR_CT1_meta_THR2","CR_CT1_meta_THR3","CR_CT1_meta_THR4" )){
  print(ggplot()+
    stat_density_2d(aes(CR_CT1_meta_back$Y,-(CR_CT1_meta_back$X)), 
                  geom = "polygon", colour="light grey",contour_var = "count")+
    stat_density_2d(aes(get(x)$Y,-(get(x)$X) , fill = stat(piece)), 
                  geom = "polygon",contour_var = "ndensity")+
    scale_fill_gradientn(colors = plasma(10),limits=c(0,20), breaks=c(0, 5, 10 ,15,20),labels=c(0, 5, 10 ,15,20))+
    labs(x= "x", y = "y")+
    xlim(-30, 180)+
    ylim(-90, 20)+
    theme_bw()+ggtitle(a[b],as.character(length(get(x)$X))))
  b=b+1
}


a <- c("CR_SP1 > THR1, nr spots:","CR_SP1 > THR2, nr spots:","CR_SP1 > THR3, nr spots:","CR_SP1 > THR4, nr spots")
b=1
for (x in c("CR_SP1_meta_THR1","CR_SP1_meta_THR2","CR_SP1_meta_THR3","CR_SP1_meta_THR4" )){
  print(ggplot()+
    stat_density_2d(aes(CR_SP1_meta_back$Y,-(CR_SP1_meta_back$X)), 
                  geom = "polygon", colour="light grey",contour_var = "count")+
    stat_density_2d(aes(get(x)$Y,-(get(x)$X) , fill = stat(piece)), 
                  geom = "polygon",contour_var = "ndensity")+
    scale_fill_gradientn(colors = plasma(10),limits=c(0,20), breaks=c(0, 5, 10 ,15,20),labels=c(0, 5, 10 ,15,20))+
    labs(x= "x", y = "y")+
    xlim(-30, 180)+
    ylim(-90, 20)+
    theme_bw()+ggtitle(a[b],as.character(length(get(x)$X)))+
    stat_density_2d(aes(CR_SP1_density_meta$Y,-(CR_SP1_density_meta$X) , fill = stat(piece)), 
                  geom = "polygon", colour="green",contour_var = "count"))
  b=b+1
}

a <- c("HE_CT1 > THR1, nr spots:","HE_CT1 > THR2, nr spots:","HE_CT1 > THR3, nr spots:","HE_CT1 > THR4, nr spots")
b=1
for (x in c("HE_CT1_meta_THR1","HE_CT1_meta_THR2","HE_CT1_meta_THR3","HE_CT1_meta_THR4" )){
  print(ggplot()+
    stat_density_2d(aes(HE_CT1_meta_back$Y,-(HE_CT1_meta_back$X)), 
                  geom = "polygon", colour="light grey",contour_var = "count")+
    stat_density_2d(aes(get(x)$Y,-(get(x)$X) , fill = stat(piece)), 
                  geom = "polygon",contour_var = "ndensity")+
    scale_fill_gradientn(colors = plasma(10),limits=c(0,20), breaks=c(0, 5, 10 ,15,20),labels=c(0, 5, 10 ,15,20))+
    labs(x= "x", y = "y")+
    xlim(-30, 180)+
    ylim(-90, 20)+
    theme_bw()+ggtitle(a[b],as.character(length(get(x)$X))))
  b=b+1
}

a <- c("HE_CT2 > THR1, nr spots:","HE_CT2 > THR2, nr spots:","HE_CT2 > THR3, nr spots:","HE_CT2 > THR4, nr spots")
b=1
for (x in c("HE_CT2_meta_THR1","HE_CT2_meta_THR2","HE_CT2_meta_THR3","HE_CT2_meta_THR4" )){
  print(ggplot()+
    stat_density_2d(aes(HE_CT2_meta_back$Y,-(HE_CT2_meta_back$X)), 
                  geom = "polygon", colour="light grey",contour_var = "count")+
    stat_density_2d(aes(get(x)$Y,-(get(x)$X) , fill = stat(piece)), 
                  geom = "polygon",contour_var = "ndensity")+
    scale_fill_gradientn(colors = plasma(10),limits=c(0,20), breaks=c(0, 5, 10 ,15,20),labels=c(0, 5, 10 ,15,20))+
    labs(x= "x", y = "y")+
    xlim(-30, 180)+
    ylim(-90, 20)+
    theme_bw()+ggtitle(a[b],as.character(length(get(x)$X))))
  b=b+1
}


a <- c("HE_SP1 > THR1, nr spots:","HE_SP1 > THR2, nr spots:")
b=1
for (x in c("HE_SP1_meta_THR1","HE_SP1_meta_THR2")){
  print(ggplot()+
    stat_density_2d(aes(HE_SP1_meta_back$Y,-(HE_SP1_meta_back$X)), 
                  geom = "polygon", colour="light grey",contour_var = "count")+
    stat_density_2d(aes(get(x)$Y,-(get(x)$X) , fill = stat(piece)), 
                  geom = "polygon",contour_var = "ndensity")+
    scale_fill_gradientn(colors = plasma(10),limits=c(0,20), breaks=c(0, 5, 10 ,15,20),labels=c(0, 5, 10 ,15,20))+
    labs(x= "x", y = "y")+
    xlim(-30, 180)+
    ylim(-90, 20)+
    theme_bw()+ggtitle(a[b],as.character(length(get(x)$X))))
  b=b+1
}


plot.new()                                  # Create empty plot


knitr::include_graphics("D:/Poland/PHD/spatial/Second_set/Untitled.jpg")  


a <- c("HE_SP2 > THR1, nr spots:","HE_SP2 > THR2, nr spots:","HE_SP2 > THR3, nr spots:","HE_SP2 > THR4, nr spots")
b=1
for (x in c("HE_SP2_meta_THR1","HE_SP2_meta_THR2","HE_SP2_meta_THR3","HE_SP2_meta_THR4" )){
  print(ggplot()+
    stat_density_2d(aes(HE_SP2_meta_back$Y,-(HE_SP2_meta_back$X)), 
                  geom = "polygon", colour="light grey",contour_var = "count")+
    stat_density_2d(aes(get(x)$Y,-(get(x)$X) , fill = stat(piece)), 
                  geom = "polygon",contour_var = "ndensity")+
    scale_fill_gradientn(colors = plasma(10),limits=c(0,20), breaks=c(0, 5, 10 ,15,20),labels=c(0, 5, 10 ,15,20))+
    labs(x= "x", y = "y")+
    xlim(-30, 180)+
    ylim(-90, 20)+
    theme_bw()+ggtitle(a[b],as.character(length(get(x)$X))))
  b=b+1
}


for (x in c("THR1", "THR2", "THR3", "THR4")){
  CR_control_11 <- subset(CR_control_1, subset = percent.AD > get(x))
print(SpatialDimPlot(CR_control_11, label = FALSE, label.size =2, crop = FALSE))
}

for (x in c("THR1", "THR2", "THR3", "THR4")){
  CR_sample_11 <- subset(CR_sample_1, subset = percent.AD > get(x))
print(SpatialDimPlot(CR_sample_11, label = FALSE, label.size =2, crop = FALSE))
}

for (x in c("THR1", "THR2", "THR3", "THR4")){
  HE_control_11 <- subset(HE_control_1, subset = percent.AD > get(x))
print(SpatialDimPlot(HE_control_11, label = FALSE, label.size =2, crop = FALSE))
}


for (x in c("THR1", "THR2", "THR3", "THR4")){
  HE_control_22 <- subset(HE_control_2, subset = percent.AD > get(x))
print(SpatialDimPlot(HE_control_22, label = FALSE, label.size =2, crop = FALSE))
}

for (x in c("THR1", "THR2")){
  HE_sample_11 <- subset(HE_sample_1, subset = percent.AD > get(x))
print(SpatialDimPlot(HE_sample_11, label = FALSE, label.size =2, crop = FALSE))
}


plot.new()                                  # Create empty plot

knitr::include_graphics("D:/Poland/PHD/spatial/Second_set/Untitled.jpg")  

for (x in c("THR1", "THR2", "THR3", "THR4")){
  HE_sample_22 <- subset(HE_sample_2, subset = percent.AD > get(x))
print(SpatialDimPlot(HE_sample_22, label = FALSE, label.size =2, crop = FALSE))
}

```
Some instances of the resluts: 
In the background, you can observe the intact image histology area and in the front, the distribution of gene expression.
![Spatial distribution of gene expressions](https://github.com/ElyasMo/Spatial_distirbution_gene_expression/blob/main/Distributions.PNG)

```{r , warning = FALSE, message = FALSE, cache=TRUE, fig.show='hold', out.width="100%", fig.fullwidth = TRUE}
boxplot(CR_sample_11@meta.data$percent.AD, CR_control_11@meta.data$percent.AD, HE_sample_11@meta.data$percent.AD, HE_control_11@meta.data$percent.AD, HE_sample_22@meta.data$percent.AD, HE_control_22@meta.data$percent.AD , main = "% of remained desired gene after filtration", names=c("CR_SP1", "CR_CT1", "HE_SP1", "HE_CT1", "HE_SP2", "HE_CT2"), col = c( "dark blue","orange","red", "yellow", "light blue", "purple"), las = 2)
```
![The box plot of the remained % of desired gene](https://github.com/ElyasMo/Spatial_distirbution_gene_expression/blob/main/Boxplot2.png)


