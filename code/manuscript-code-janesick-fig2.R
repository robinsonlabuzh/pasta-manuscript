library(dplyr)
library(scran)
library(spdep)
library(tidyr)
library(ggplot2)
library(Voyager)
library(SFEData)
library(spatstat)
library(openxlsx)
library(spatialFDA)
library(BiocParallel)
library(SpatialExperiment)
library(SpatialFeatureExperiment)
library(patchwork)
library(pals)

spe <- STexampleData::Janesick_breastCancer_Xenium_rep1()

sfe <- toSpatialFeatureExperiment(spe)

### code adapted from https://lmweber.org/OSTA/pages/crs-spat-stat.html ###
### written by Samuel Gunz and Martin Emons 2025 ###

# load the official 10X annotations
labels <- read.xlsx("https://cdn.10xgenomics.com/raw/upload/v1695234604/Xenium%20Preview%20Data/Cell_Barcode_Type_Matrices.xlsx", sheet = 4)

labels$cell_id <- (labels$Barcode)

# add the cell type labels to the spe
matchedDf <- as.data.frame(colData(sfe)) |>
  left_join(as.data.frame(labels), by = join_by("cell_id"))

colData(sfe) <- DataFrame(matchedDf)

sfe <- sfe[, colSums(counts(sfe)) > 0]
sfe <- sfe[, !is.na(sfe$Cluster)]
sfe <- logNormCounts(sfe)

xy <- spatialCoords(sfe)

colGraph(sfe, "knn6") <-
  findSpatialNeighbors(
    sfe,
    type = "centroids",
    method = "knearneigh",
    k = 6
  )

plotColGraph(sfe,
             colGraphName = "knn6",
             colGeometryName = "centroids"
) + theme_void()

# get all gene probes
geneProbes <- rowData(sfe)[rowData(sfe)$Type == "Gene Expression", "Symbol"]
# select receptor genes
receptorGenes <- c("ERBB2", "ESR1", "PGR")

sfe <- runUnivariate(sfe,
                     type = "moran",
                     features = receptorGenes,
                     colGraphName = "knn6"
)

p0 <- plotSpatialFeature(sfe, receptorGenes, ncol = 3,
                         size = 1E-5, alpha = 0.75, scattermore = FALSE)



sfe <- runUnivariate(sfe,
                     type = "localmoran",
                     features = receptorGenes,
                     colGraphName = "knn6"
)

p <- plotLocalResult(sfe,
                     name = "localmoran",
                     features = receptorGenes,
                     colGeometryName = "centroids",
                     divergent = TRUE,
                     diverge_center = 0,
                     size = 0.001,
                     ncol = 3,
                     alpha = 0.75,
                     scattermore = FALSE)

hvgs <- getTopHVGs(sfe, n = 20)

res <- calculateBivariate(sfe,
                          type = "lee",
                          feature1 = receptorGenes,
                          colGraphName =  "knn6"
)

sfe <- runBivariate(
  sfe,
  "locallee",
  colGraphName =  "knn6",
  feature1 = receptorGenes
)

p1 <- plotLocalResult(
  sfe,
  name = "locallee",
  features = "ERBB2__ESR1",
  colGeometryName = "centroids",
  divergent = TRUE,
  diverge_center = 0,
  size = 0.001,
  alpha = 0.75,
  scattermore = FALSE
) + ggtitle(NULL)

p2 <- plotLocalResult(
  sfe,
  name = "locallee",
  features = "ERBB2__PGR",
  colGeometryName = "centroids",
  divergent = TRUE,
  diverge_center = 0,
  size = 0.001,
  alpha = 0.75,
  scattermore = FALSE
) + ggtitle(NULL)

p3 <- plotLocalResult(
  sfe,
  name = "locallee",
  features = "ESR1__PGR",
  colGeometryName = "centroids",
  divergent = TRUE,
  diverge_center = 0,
  size = 0.001,
  alpha = 0.75,
  scattermore = FALSE
) + ggtitle(NULL)


p0 <- p0 + geom_point(size = 1E-5)
p1 <- p1 + geom_point(size = 1E-5)
p2 <- p2 + geom_point(size = 1E-5)
p3 <- p3 + geom_point(size = 1E-5)

pLee <- wrap_plots(list(p1,p2,p3))

sfe <- runMultivariate(sfe, type = "localC_perm_multi",
                       subset_row = receptorGenes,
                       colGraphName = "knn6")

colData(sfe)$logLocalC <- log(reducedDim(sfe, "localC_perm_multi")["localC_perm_multi"])


reducedDim(sfe, "localC_perm_multi")["log_localC_perm_multi"] <- log(reducedDim(sfe, "localC_perm_multi")["localC_perm_multi"])

# stored as spatially reduced dim; plot it in this way
p4 <- spatialReducedDim(sfe, "localC_perm_multi", c(13, 11),
                        size = 0.001, alpha = 0.75, scattermore = FALSE, divergent = TRUE,
                        diverge_center = 0)
p4 <- p4 + geom_point(size = 1E-5)

p4

pAll <- p0/p/pLee /p4
pAll <- pAll + plot_annotation(tag_levels = 'A')

#pAll

ggsave(plot =pAll, "outs/fig2.png", width = 10, height = 10, dpi = 150)


#save triple positive result - multivar Geary's C
p_ERBB2ESR1PGR <- p4[[1]] + theme(legend.title = element_blank()) +
  scale_color_gradientn(colours = brewer.piyg(100), guide = "colourbar") +
  ggtitle('Spatial Correlation ERBB2, ESR1, PGR') +
  geom_density_2d()
ggsave(plot = p_ERBB2ESR1PGR, "outs/fig3+.pdf", width = 5, height = 5)


#save double positive result - bivar Lee's L
p_ERBB2ESR1 <- p1  + 
  scale_color_gradientn(colours = brewer.rdbu(100), guide = "colourbar") +
  theme(legend.title = element_blank())+
  ggtitle('Spatial Correlation ERBB2, ESR1')
ggsave(plot =p_ERBB2ESR1, "outs/fig2+.pdf", width = 5, height = 5)

#save single positive result - local Moran's I
p_ERBB2 <- p[[1]]  + 
  scale_color_gradientn(colours = brewer.puor(100), guide = "colourbar") +
  #scale_color_gradientn(colours = brewer.spectral(100), guide = "colourbar") +
  theme(legend.title = element_blank())+
  ggtitle('Spatial Correlation ERBB2')
ggsave(plot =p_ERBB2, "outs/fig1+.pdf", width = 5, height = 5)
