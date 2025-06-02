library(dplyr)
library(spdep)
library(tidyr)
library(ggplot2)
library(Voyager)
library(spatstat)
library(openxlsx)
library(spatialFDA)
library(BiocParallel)
library(SpatialExperiment)
library(SpatialFeatureExperiment)
library(pals)
library(SpaNorm)
library(patchwork)

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

# perform a spatially aware normalisation
sfe <- SpaNorm(sfe, sample.p = 0.01, tol = 1e-03)

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


### single positive analysis ###

p0 <- plotSpatialFeature(sfe, receptorGenes, ncol = 1,
                         size = 1E-5, alpha = 0.75, scattermore = FALSE)

sfe <- runUnivariate(sfe,
                     type = "localmoran_perm",
                     features = receptorGenes,
                     colGraphName = "knn6"
)

#change levels
localResult(sfe, feature = receptorGenes)$mean <-
  factor(localResult(sfe, feature = receptorGenes)$mean, 
         level = c( "High-High", "Low-Low" ,  "High-Low",  "Low-High"))

pUnivar <- plotLocalResult(sfe,
                     name = "localmoran_perm",
                     features = receptorGenes,
                     attribute = "mean",
                     colGeometryName = "centroids",
                     divergent = TRUE,
                     diverge_center = 0,
                     size = 0.001,
                     ncol = 1,
                     alpha = 0.75,
                     scattermore = FALSE) 

pUnivar[[1]]  <- pUnivar[[1]] + guides(color = guide_legend(override.aes = list(size = 3)))
pUnivar[[2]]  <- pUnivar[[2]] + guides(color = guide_legend(override.aes = list(size = 3)))
pUnivar[[3]]  <- pUnivar[[3]] + guides(color = guide_legend(override.aes = list(size = 3)))

### double positive analysis ###

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

pBivar <- wrap_plots(list(p1,p2,p3), nrow = 3)


### triple positive analysis ###

sfe <- runMultivariate(sfe, type = "localC_perm_multi",
                       nsim = 499, 
                       subset_row = receptorGenes,
                       colGraphName = "knn6",
                       exprs_values = "logcounts")

#create a 4x4 clustering table as in https://mkram01.github.io/EPI563-SpatialEPI/spatial-structure-and-clustering-i-morans-i-and-lisa.html
#code adapted from the webpage lincensed under CC-BY-NC-SA

reducedDim(sfe, "localC_perm_multi") <- reducedDim(sfe, "localC_perm_multi") %>% 
  mutate(
    auto_cluster = factor(case_when(
      cluster == "Positive" & `-log10p_adj Sim` >= -log10(0.05) ~ 'Positive',
      cluster == "Negative" & `-log10p_adj Sim` >= -log10(0.05) ~ 'Negative',
      `-log10p_adj Sim` < -log10(0.05) ~ 'Non-significant',),
      levels = c('Positive', 'Negative', 'Non-significant')),
    manual_cluster = factor(case_when(
      localC_perm_multi < 1 & `-log10p_adj Sim` >= -log10(0.05) ~ 'Positive',
      localC_perm_multi >= 1 & `-log10p_adj Sim` >= -log10(0.05) ~ 'Negative',
      `-log10p_adj Sim` < -log10(0.05) ~ 'Non-significant',),
      levels = c('Positive', 'Negative', 'Non-significant'))
  )


### end ###

# stored as spatially reduced dim; plot it in this way
pMultivar <- spatialReducedDim(sfe, "localC_perm_multi", c(1, 11, 12),
                        size = 0.001, alpha = 0.75, scattermore = FALSE, divergent = TRUE,
                        diverge_center = 0, ncol = 1)
pMultivar <- pMultivar + geom_point(size = 1E-5) + guides(color = guide_legend(override.aes = list(size = 3)))

pMultivar

pDens <- cbind(reducedDim(sfe, "localC_perm_multi"), spatialCoords(sfe)) |>
  filter(manual_cluster == "Positive") |>
  ggplot(aes(x = x_centroid, y = y_centroid)) +
  geom_point(size = 0.01) +
  theme_void() +
  geom_density_2d(alpha = 0.5) +
  labs(title = "Density of positive & significant points") +
  coord_equal()

pMultivar[[3]] + scale_color_manual(values = c("blue","red", "grey")) + 
  guides(color = guide_legend(override.aes = list(size = 3))) +
  pDens

pAlla <- p0|pUnivar
pAlla <- pAlla +  plot_annotation(tag_levels = 'A')

pAllb <- pBivar|pMultivar
pAllb <- pAllb +  plot_annotation(tag_levels = 'A')

ggsave(plot =pAlla, "outs/fig2a.png", width = 10, height = 10, dpi = 200)
ggsave(plot =pAllb, "outs/fig2b.png", width = 10, height = 10, dpi = 200)



### consensus analysis on the single positive Geary's C regions ###

# alternative approach by finding the consensus of the single positive regions
df <- data.frame()

df <- lapply(receptorGenes, function(elem){
  localResult(sfe, feature = elem)
}) %>% as.data.frame() %>%
  cbind()


df <- df %>%
  mutate(consensus = factor(case_when(
    mean == "High-High" & mean.1 == "High-High" & mean.2 == "High-High" ~ "ERBB2+/ESR1+/PGR+",
    mean == "High-High" & mean.1 == "High-High" ~ "ERBB2+/ESR1+",
    mean == "High-High" & mean.2 == "High-High" ~ "ERBB2+/PGR+",
    mean.1 == "High-High" & mean.2 == "High-High" ~ "ESR1+/PGR+",
    mean == "High-High" ~ "ERBB2+",
    mean.1 == "High-High" ~ "ESR1+",
    mean.2 == "High-High" ~ "PGR+",
    X.log10p_adj.Sim < -log(0.05) ~ "Non-significant",
    X.log10p_adj.Sim.1 < -log(0.05) ~ "Non-significant",
    X.log10p_adj.Sim.2 < -log(0.05) ~ "Non-significant",
    .default = "Not applicable",
  ), levels = c("ERBB2+/ESR1+/PGR+", "ERBB2+/ESR1+", "ERBB2+/PGR+", "ESR1+/PGR+",
                "ERBB2+", "ESR1+", "PGR+", "Non-significant", "Not applicable")
  ))

df <- df %>% cbind(xy)
pConsensus <- ggplot(df, aes(x_centroid, y_centroid, color = consensus)) +
  geom_point(size = 1E-8) +
  theme_light() + 
  guides(color = guide_legend(override.aes = list(size = 3), 
                              theme = theme(
                                legend.text = element_text(size = 15),
                                legend.title = element_text(size = 15)))) +
  coord_equal() +
  scale_color_brewer(palette = "Accent") +
  ggtitle("Consensus of Moran's Scatterplot for *ERBB2*, *ESR1*, *PGR*")+
  theme(plot.title = ggtext::element_markdown(size = 16))

pConsensus <- ggrastr::rasterize(pConsensus, layers='Point', dpi=300)

ggsave(plot = pConsensus, "outs/consensus.pdf", width = 15, height = 10)

