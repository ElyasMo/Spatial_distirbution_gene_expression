# Spatial_distirbution_gene_expression
Here I have developed a pipeline to monitor the spatial distribution of gene expression in various Seurat objects based on four different thresholds. 

**The R markdown output**
---
title: "ITGAM(CD11b)-spatial distribution of expression"
author: "Elyas Mohammadi"
date: "April 2, 2021"
output: html_document
---

```{r, warning = FALSE, message = FALSE, cache=TRUE}

Genes <- ("ITGAM")
THR1 = 0
THR2 = 0.02
THR3 = 0.04
THR4 = 0.05


library(Seurat)
library(sctransform)

#####Reading in the CR_sample_2
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
setwd("D:/Poland/PHD/spatial/project/CR/CR_control_2339/")
Directory="D:/Poland/PHD/spatial/project/CR/CR_control_2339/"
CR_control_1<-Load10X_Spatial(Directory,filename = "filtered_feature_bc_matrix.h5",
                              assay = "Spatial",
                              slice = "CRctr1",
                              filter.matrix = TRUE,
                              to.upper = FALSE)



##Reading the 10xgenomics AD genes
#Genes <- read.csv("D:/Poland/PHD/spatial/Second_set/DEGs/AD_genes_10xgenomics.csv", row.names=NULL)
#Genes <- as.list(Genes[,1])


##Calculating the percentage of AD genes in each spot
AD_gene <- match(Genes,rownames(CR_sample_1@assays$Spatial))
AD_genes <- rownames(CR_sample_1@assays$Spatial)[AD_gene]
AD_genes <- AD_genes[!is.na(AD_genes)]
CR_sample_1[["percent.AD"]]<-PercentageFeatureSet(CR_sample_1,features=AD_genes)

AD_gene <- match(Genes,rownames(CR_control_1@assays$Spatial))
AD_genes <- rownames(CR_control_1@assays$Spatial)[AD_gene]
AD_genes <- AD_genes[!is.na(AD_genes)]
CR_control_1[["percent.AD"]]<-PercentageFeatureSet(CR_control_1,features=AD_genes)

AD_gene <- match(Genes,rownames(HE_sample_1@assays$Spatial))
AD_genes <- rownames(HE_sample_1@assays$Spatial)[AD_gene]
AD_genes <- AD_genes[!is.na(AD_genes)]
HE_sample_1[["percent.AD"]]<-PercentageFeatureSet(HE_sample_1,features=AD_genes)

AD_gene <- match(Genes,rownames(HE_control_1@assays$Spatial))
AD_genes <- rownames(HE_control_1@assays$Spatial)[AD_gene]
AD_genes <- AD_genes[!is.na(AD_genes)]
HE_control_1[["percent.AD"]]<-PercentageFeatureSet(HE_control_1,features=AD_genes)

AD_gene <- match(Genes,rownames(HE_sample_2@assays$Spatial))
AD_genes <- rownames(HE_sample_2@assays$Spatial)[AD_gene]
AD_genes <- AD_genes[!is.na(AD_genes)]
HE_sample_2[["percent.AD"]]<-PercentageFeatureSet(HE_sample_2,features=AD_genes)

AD_gene <- match(Genes,rownames(HE_control_2@assays$Spatial))
AD_genes <- rownames(HE_control_2@assays$Spatial)[AD_gene]
AD_genes <- AD_genes[!is.na(AD_genes)]
HE_control_2[["percent.AD"]]<-PercentageFeatureSet(HE_control_2,features=AD_genes)


##adding the percentage of MT genes to metadata
CR_control_1[["percent.mt"]] <- PercentageFeatureSet(CR_control_1, pattern = "^MT-")
CR_sample_1[["percent.mt"]] <- PercentageFeatureSet(CR_sample_1, pattern = "^MT-")
HE_control_1[["percent.mt"]] <- PercentageFeatureSet(HE_control_1, pattern = "^MT-")
HE_sample_1[["percent.mt"]] <- PercentageFeatureSet(HE_sample_1, pattern = "^MT-")
HE_control_2[["percent.mt"]] <- PercentageFeatureSet(HE_control_2, pattern = "^MT-")
HE_sample_2[["percent.mt"]] <- PercentageFeatureSet(HE_sample_2, pattern = "^MT-")



##Filtering the metadata basd on various criteria including percent.AD > 1.5
CR_control_1 <- subset(CR_control_1, subset = nFeature_Spatial > 200 & nFeature_Spatial < 7000 & percent.mt < 15)
CR_sample_1 <- subset(CR_sample_1, subset = nFeature_Spatial > 200 & nFeature_Spatial < 7000 & percent.mt < 15 )
HE_control_1 <- subset(HE_control_1, subset = nFeature_Spatial > 200 & nFeature_Spatial < 7000 & percent.mt < 15)
HE_sample_1 <- subset(HE_sample_1, subset = nFeature_Spatial > 200 & nFeature_Spatial < 7000 & percent.mt < 15)
HE_control_2 <- subset(HE_control_2, subset = nFeature_Spatial > 200 & nFeature_Spatial < 7000 & percent.mt < 15)
HE_sample_2 <- subset(HE_sample_2, subset = nFeature_Spatial > 200 & nFeature_Spatial < 7000 & percent.mt < 15)


##A good estimation for the threasholds
quantile(CR_sample_1@meta.data$percent.AD,prob=seq(0,1,.1))
quantile(CR_control_1@meta.data$percent.AD,prob=seq(0,1,.1))
quantile(HE_control_1@meta.data$percent.AD,prob=seq(0,1,.1))
quantile(HE_sample_1@meta.data$percent.AD,prob=seq(0,1,.1))
quantile(HE_control_2@meta.data$percent.AD,prob=seq(0,1,.1))
quantile(HE_sample_2@meta.data$percent.AD,prob=seq(0,1,.1))

```

**Adding the cordinates**
The metadata is filtered and we would like to see the locations which have the accumulated spots with high percentage of determined genes

```{r , warning = FALSE, message = FALSE, cache=TRUE}
##inserting the cordinates
setwd("D:/Poland/PHD/spatial/Second_set/cordinations/")
CR_CT1_spatial <- read.csv("D:/Poland/PHD/spatial/Second_set/cordinations/CR_CT1_spatial.csv", row.names=1)
CR_SP1_spatial <- read.csv("D:/Poland/PHD/spatial/Second_set/cordinations/CR_SP1_spatial..csv", row.names=1)
HE_CT1_spatial <- read.csv("D:/Poland/PHD/spatial/Second_set/cordinations/HE_CT1_spatial.csv", row.names=1)
HE_CT2_spatial <- read.csv("D:/Poland/PHD/spatial/Second_set/cordinations/HE_CT2_spatial.csv", row.names=1)
HE_SP1_spatial <- read.csv("D:/Poland/PHD/spatial/Second_set/cordinations/HE_SP1_spatial.csv", row.names=1)
HE_SP2_spatial <- read.csv("D:/Poland/PHD/spatial/Second_set/cordinations/HE_SP2_spatial.csv", row.names=1)

##Adding cordinates to the metadata
CR_CT1_meta <- cbind(CR_control_1@meta.data[intersect(rownames(CR_control_1@meta.data), rownames(CR_CT1_spatial)),], CR_CT1_spatial[intersect(rownames(CR_control_1@meta.data), rownames(CR_CT1_spatial)),])

CR_SP1_meta <-cbind(CR_sample_1@meta.data[intersect(rownames(CR_sample_1@meta.data), rownames(CR_SP1_spatial)),], CR_SP1_spatial[intersect(rownames(CR_sample_1@meta.data), rownames(CR_SP1_spatial)),])

HE_CT1_meta <-cbind(HE_control_1@meta.data[intersect(rownames(HE_control_1@meta.data), rownames(HE_CT1_spatial)),], HE_CT1_spatial[intersect(rownames(HE_control_1@meta.data), rownames(HE_CT1_spatial)),])

HE_CT2_meta <- cbind(HE_control_2@meta.data[intersect(rownames(HE_control_2@meta.data), rownames(HE_CT2_spatial)),], HE_CT2_spatial[intersect(rownames(HE_control_2@meta.data), rownames(HE_CT2_spatial)),])

HE_SP1_meta <- cbind(HE_sample_1@meta.data[intersect(rownames(HE_sample_1@meta.data), rownames(HE_SP1_spatial)),], HE_SP1_spatial[intersect(rownames(HE_sample_1@meta.data), rownames(HE_SP1_spatial)),])

HE_SP2_meta <- cbind(HE_sample_2@meta.data[intersect(rownames(HE_sample_2@meta.data), rownames(HE_SP2_spatial)),],  HE_SP2_spatial[intersect(rownames(HE_sample_2@meta.data), rownames(HE_SP2_spatial)),])

```
**Modifying the metadata for manual grouping in order to illustrate the amyloid spots**
```{r, warning = FALSE, message = FALSE, cache=TRUE}
CR_sample_1 <- SCTransform(CR_sample_1, assay = "Spatial", verbose = FALSE)
CR_sample_1<- FindVariableFeatures(CR_sample_1, selection.method = "vst", 
                                   nfeatures = 3000, verbose = FALSE)


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

SpatialDimPlot(CR_sample_1, label = FALSE, label.size =2, crop = FALSE)


```

**Plotting**

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

**In this section you can see CR_control_1 plots based on the threashold you determined at the beginning (THR threasholds)**

```{r , warning = FALSE, message = FALSE, cache=TRUE, fig.width=6, fig.show='hold', out.width="50%"}

ggplot()+
  stat_density_2d(aes(CR_CT1_meta_THR1$Y,-(CR_CT1_meta_THR1$X) , fill = after_stat(level)), 
                  geom = "polygon")+
  scale_fill_gradientn(colors = plasma(10))+
  labs(x= "x", y = "y", title = "2D KDE plot for CR_CT1 with determined gene > THR1")+
  xlim(-10, 140)+
  ylim(-80, -1)+
  theme_bw()

CR_control_1 <- subset(CR_control_1, subset = percent.AD > THR1)
SpatialDimPlot(CR_control_1, label = FALSE, label.size =2, crop = FALSE)



ggplot()+
  stat_density_2d(aes(CR_CT1_meta_THR2$Y,-(CR_CT1_meta_THR2$X) , fill = after_stat(level)), 
                  geom = "polygon")+
  scale_fill_gradientn(colors = plasma(10))+
  labs(x= "x", y = "y", title = "2D KDE plot for CR_CT1 with determined genes > THR2")+
  xlim(-10, 180)+
  ylim(-90, 5)+
  theme_bw()

CR_control_1 <- subset(CR_control_1, subset = percent.AD > THR2)
SpatialDimPlot(CR_control_1, label = FALSE, label.size =2, crop = FALSE)



ggplot()+
  stat_density_2d(aes(CR_CT1_meta_THR3$Y,-(CR_CT1_meta_THR3$X) , fill = after_stat(level)), 
                  geom = "polygon")+
  scale_fill_gradientn(colors = plasma(10))+
  labs(x= "x", y = "y", title = "2D KDE plot for CR_CT1 with determined genes > THR3")+
  xlim(-10, 170)+
  ylim(-90, 5)+
  theme_bw()

CR_control_1 <- subset(CR_control_1, subset = percent.AD > THR3)
SpatialDimPlot(CR_control_1, label = FALSE, label.size =2, crop = FALSE)



ggplot()+
  stat_density_2d(aes(CR_CT1_meta_THR4$Y,-(CR_CT1_meta_THR4$X) , fill = after_stat(level)), 
                  geom = "polygon")+
  scale_fill_gradientn(colors = plasma(10))+
  labs(x= "x", y = "y", title = "2D KDE plot for CR_CT1 with determined genes > THR4")+
  xlim(-10, 170)+
  ylim(-90, 20)+
  theme_bw()

CR_control_1 <- subset(CR_control_1, subset = percent.AD > THR4)
SpatialDimPlot(CR_control_1, label = FALSE, label.size =2,crop = FALSE)


```


**In this section you can see CR_sample_1 plots based on the threashold you determined at the beginning (THR threasholds)**
```{r , warning = FALSE, message = FALSE, cache=TRUE, fig.width=6, fig.show='hold', out.width="50%"}
ggplot()+
  stat_density_2d(aes(CR_SP1_meta_THR1$Y,-(CR_SP1_meta_THR1$X) , fill = after_stat(level)), 
                  geom = "polygon")+
  scale_fill_gradientn(colors = plasma(10))+
  labs(x= "x", y = "y", title = "2D KDE plot for CR_SP1 with determined gene > THR1")+
  xlim(-10, 140)+
  ylim(-80, -1)+
  theme_bw()

CR_sample_1 <- subset(CR_sample_1, subset = percent.AD > THR1)
SpatialDimPlot(CR_sample_1, label = FALSE, label.size =2,crop = FALSE)


ggplot()+
  stat_density_2d(aes(CR_SP1_meta_THR2$Y,-(CR_SP1_meta_THR2$X) , fill = after_stat(level)), 
                  geom = "polygon")+
  scale_fill_gradientn(colors = plasma(10))+
  labs(x= "x", y = "y", title = "2D KDE plot for CR_SP1 with determined genes > THR2")+
  xlim(-10, 130)+
  ylim(-75, 20)+
  theme_bw()

CR_sample_1 <- subset(CR_sample_1, subset = percent.AD > THR2)

SpatialDimPlot(CR_sample_1, label = FALSE, label.size =2, crop = FALSE)


ggplot()+
  stat_density_2d(aes(CR_SP1_meta_THR3$Y,-(CR_SP1_meta_THR3$X) , fill = after_stat(level)), 
                  geom = "polygon")+
  scale_fill_gradientn(colors = plasma(10))+
  labs(x= "x", y = "y", title = "2D KDE plot for CR_SP1 with determined genes > THR3")+
  xlim(-20, 140)+
  ylim(-75, 20)+
  theme_bw()

CR_sample_1 <- subset(CR_sample_1, subset = percent.AD > THR3)

SpatialDimPlot(CR_sample_1, label = FALSE, label.size =2, crop = FALSE)





ggplot()+
  stat_density_2d(aes(CR_SP1_meta_THR4$Y,-(CR_SP1_meta_THR4$X) , fill = after_stat(level)), 
                  geom = "polygon")+
  scale_fill_gradientn(colors = plasma(10))+
  labs(x= "x", y = "y", title = "2D KDE plot for CR_SP1 with determined genes > THR4")+
  xlim(20, 140)+
  ylim(-40, 15)+
  theme_bw()

CR_sample_1 <- subset(CR_sample_1, subset = percent.AD > THR4)
SpatialDimPlot(CR_sample_1, label = FALSE, label.size =2, crop = FALSE)
```
**In this section you can see HE_control_1 plots based on the threashold you determined at the beginning (THR threasholds)**

```{r , warning = FALSE, message = FALSE, cache=TRUE, fig.width=6, fig.show='hold', out.width="50%"}

ggplot()+
  stat_density_2d(aes(HE_CT1_meta_THR1$Y,-(HE_CT1_meta_THR1$X) , fill = after_stat(level)), 
                  geom = "polygon")+
  scale_fill_gradientn(colors = plasma(10))+
  labs(x= "x", y = "y", title = "2D KDE plot for HE_CT1 with determined gene > THR1")+
  xlim(-30, 170)+
  ylim(-80, 15)+
  theme_bw()

HE_control_1 <- subset(HE_control_1, subset = percent.AD > THR1)
SpatialDimPlot(HE_control_1, label = FALSE, label.size =2, crop = FALSE)



ggplot()+
  stat_density_2d(aes(HE_CT1_meta_THR2$Y,-(HE_CT1_meta_THR2$X) , fill = after_stat(level)), 
                  geom = "polygon")+
  scale_fill_gradientn(colors = plasma(10))+
  labs(x= "x", y = "y", title = "2D KDE plot for HE_CT1 with determined genes > THR2")+
  xlim(-30, 130)+
  ylim(-90, 15)+
  theme_bw()

HE_control_1 <- subset(HE_control_1, subset = percent.AD > THR2)
SpatialDimPlot(HE_control_1, label = FALSE, label.size =2, crop = FALSE)



ggplot()+
  stat_density_2d(aes(HE_CT1_meta_THR3$Y,-(HE_CT1_meta_THR3$X) , fill = after_stat(level)), 
                  geom = "polygon")+
  scale_fill_gradientn(colors = plasma(10))+
  labs(x= "x", y = "y", title = "2D KDE plot for HE_CT1 with determined genes > THR3")+
  xlim(-40, 130)+
  ylim(-90, 20)+
  theme_bw()

HE_control_1 <- subset(HE_control_1, subset = percent.AD > THR3)
SpatialDimPlot(HE_control_1, label = FALSE, label.size =2, crop = FALSE)



ggplot()+
  stat_density_2d(aes(HE_CT1_meta_THR4$Y,-(HE_CT1_meta_THR4$X) , fill = after_stat(level)), 
                  geom = "polygon")+
  scale_fill_gradientn(colors = plasma(10))+
  labs(x= "x", y = "y", title = "2D KDE plot for HE_CT1 with determined genes > THR4")+
  xlim(-30, 130)+
  ylim(-90, 20)+
  theme_bw()

HE_control_1 <- subset(HE_control_1, subset = percent.AD > THR4)
SpatialDimPlot(HE_control_1, label = FALSE, label.size =2,crop = FALSE)


```
**In this section you can see HE_control_2 plots based on the threashold you determined at the beginning (THR threasholds)**

```{r , warning = FALSE, message = FALSE, cache=TRUE, fig.width=6, fig.show='hold', out.width="50%"}

ggplot()+
  stat_density_2d(aes(HE_CT2_meta_THR1$Y,-(HE_CT2_meta_THR1$X) , fill = after_stat(level)), 
                  geom = "polygon")+
  scale_fill_gradientn(colors = plasma(10))+
  labs(x= "x", y = "y", title = "2D KDE plot for HE_CT1 with determined gene > THR1")+
  xlim(-30, 170)+
  ylim(-80, 15)+
  theme_bw()

HE_control_2 <- subset(HE_control_2, subset = percent.AD > THR1)
SpatialDimPlot(HE_control_2, label = FALSE, label.size =2, crop = FALSE)



ggplot()+
  stat_density_2d(aes(HE_CT2_meta_THR2$Y,-(HE_CT2_meta_THR2$X) , fill = after_stat(level)), 
                  geom = "polygon")+
  scale_fill_gradientn(colors = plasma(10))+
  labs(x= "x", y = "y", title = "2D KDE plot for HE_CT1 with determined genes > THR2")+
  xlim(-30, 130)+
  ylim(-90, 15)+
  theme_bw()

HE_control_2 <- subset(HE_control_2, subset = percent.AD > THR2)
SpatialDimPlot(HE_control_2, label = FALSE, label.size =2, crop = FALSE)



ggplot()+
  stat_density_2d(aes(HE_CT2_meta_THR3$Y,-(HE_CT2_meta_THR3$X) , fill = after_stat(level)), 
                  geom = "polygon")+
  scale_fill_gradientn(colors = plasma(10))+
  labs(x= "x", y = "y", title = "2D KDE plot for HE_CT1 with determined genes > THR3")+
  xlim(-40, 130)+
  ylim(-90, 20)+
  theme_bw()

HE_control_2 <- subset(HE_control_2, subset = percent.AD > THR3)
SpatialDimPlot(HE_control_2, label = FALSE, label.size =2, crop = FALSE)



ggplot()+
  stat_density_2d(aes(HE_CT2_meta_THR4$Y,-(HE_CT2_meta_THR4$X) , fill = after_stat(level)), 
                  geom = "polygon")+
  scale_fill_gradientn(colors = plasma(10))+
  labs(x= "x", y = "y", title = "2D KDE plot for HE_CT1 with determined genes > THR4")+
  xlim(-30, 130)+
  ylim(-90, 20)+
  theme_bw()

HE_control_2 <- subset(HE_control_2, subset = percent.AD > THR4)
SpatialDimPlot(HE_control_2, label = FALSE, label.size =2,crop = FALSE)

```
**In this section you can see HE_samlple_2 plots based on the threashold you determined at the beginning (THR threasholds)**

```{r , warning = FALSE, message = FALSE, cache=TRUE, fig.width=6, fig.show='hold', out.width="50%"}
ggplot()+
  stat_density_2d(aes(HE_SP2_meta_THR1$Y,-(HE_SP2_meta_THR1$X) , fill = after_stat(level)), 
                  geom = "polygon")+
  scale_fill_gradientn(colors = plasma(10))+
  labs(x= "x", y = "y", title = "2D KDE plot for HE_SP2 with determined gene > THR1")+
  xlim(-10, 140)+
  ylim(-90, 10)+
  theme_bw()

HE_sample_2 <- subset(HE_sample_2, subset = percent.AD > THR1)
SpatialDimPlot(HE_sample_2, label = FALSE, label.size =2,crop = FALSE)


ggplot()+
  stat_density_2d(aes(HE_SP2_meta_THR2$Y,-(HE_SP2_meta_THR2$X) , fill = after_stat(level)), 
                  geom = "polygon")+
  scale_fill_gradientn(colors = plasma(10))+
  labs(x= "x", y = "y", title = "2D KDE plot for HE_SP2 with determined genes > THR2")+
  xlim(-10, 180)+
  ylim(-90, 5)+
  theme_bw()

HE_sample_2 <- subset(HE_sample_2, subset = percent.AD > THR2)

SpatialDimPlot(HE_sample_2, label = FALSE, label.size =2, crop = FALSE)


ggplot()+
  stat_density_2d(aes(HE_SP2_meta_THR3$Y,-(HE_SP2_meta_THR3$X) , fill = after_stat(level)), 
                  geom = "polygon")+
  scale_fill_gradientn(colors = plasma(10))+
  labs(x= "x", y = "y", title = "2D KDE plot for HE_SP2 with determined genes > THR3")+
  xlim(-10, 170)+
  ylim(-90, 5)+
  theme_bw()

HE_sample_2 <- subset(HE_sample_2, subset = percent.AD > THR3)

SpatialDimPlot(HE_sample_2, label = FALSE, label.size =2, crop = FALSE)





ggplot()+
  stat_density_2d(aes(HE_SP2_meta_THR4$Y,-(HE_SP2_meta_THR4$X) , fill = after_stat(level)), 
                  geom = "polygon")+
  scale_fill_gradientn(colors = plasma(10))+
  labs(x= "x", y = "y", title = "2D KDE plot for HE_SP2 with determined genes > THR4")+
  xlim(-10, 170)+
  ylim(-90, 20)+
  theme_bw()

HE_sample_2 <- subset(HE_sample_2, subset = percent.AD > THR4)
SpatialDimPlot(HE_sample_2, label = FALSE, label.size =2, crop = FALSE)
```
