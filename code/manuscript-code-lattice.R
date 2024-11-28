library(STexampleData)
library(ggspavis)
library(viridis)
library(RColorBrewer)
library(spdep)
library(scater)
library(tmap)
library(patchwork)
library(dplyr)
library(tidyr)
library(scales)

source("code/utils.R")

spe <- Visium_mouseCoronal()
rownames(spe) <- rowData(spe)$gene_name

####
#### Source: https://lmweber.org/BestPracticesST/chapters/workflow-Visium-mouseCoronal.html
# subset to keep only spots over tissue
spe <- spe[, colData(spe)$in_tissue == 1]
dim(spe)

# identify mitochondrial genes
is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
table(is_mito)

# calculate per-spot QC metrics and store in colData
spe <- addPerCellQC(spe, subsets = list(mito = is_mito))
head(colData(spe))

hist(colData(spe)$sum, breaks = 25)
qc_lib_size <- colData(spe)$sum < 600

# plot library size vs. number of cells per spot
hist(colData(spe)$detected, breaks = 20)
qc_detected <- colData(spe)$detected < 2000


# histogram of mitochondrial read proportions
hist(colData(spe)$subsets_mito_percent, breaks = 20)

# select QC threshold for mitochondrial read proportion
qc_mito <- colData(spe)$subsets_mito_percent > 28
table(qc_mito)


discard <- qc_lib_size | qc_detected | qc_mito

colData(spe)$discard <- discard

# remove combined set of low-quality spots
spe <- spe[, !colData(spe)$discard]
dim(spe)

spe <- logNormCounts(spe)

# plot
plotVisium(spe, fill = "Nrgn", assay = "logcounts", image = FALSE)


# add some values in 'colData' to annotate spots
colData(spe)$sum <- colSums(counts(spe))
colData(spe)$logsum <- log(colSums(counts(spe)))

#saveRDS(spe, "data/Visium_mouseCorona.rds")
#create sf object

coords <- as.data.frame(spatialCoords(spe))
coords$pxl_row_in_fullres <- -coords$pxl_row_in_fullres
spsf <- st_as_sf(coords,
                 coords = c("pxl_col_in_fullres", "pxl_row_in_fullres")
)

# Select gene "Nrgn"
Nrgn <- logcounts(spe)["Nrgn", ] |>
  as.matrix() |>
  as.data.frame()

colnames(Nrgn) <- "Nrgn"

# create object
spsf <- cbind(spsf, colData(spe), Nrgn)

st_buffer(st_geometry(spsf), dist = 80) |>
  poly2nb(snap = 10) |>
  nb2listw() -> direct_neigbours

# Moran's I
moran.test(spsf[["Nrgn"]], listw = direct_neigbours)


# local moran's I
loc <- localmoran(spsf[["Nrgn"]], listw = direct_neigbours)
# extract the effect size
locEffect <- loc[, 1]
# extract the p-value and adjust for multiple testing
p.val.adj <- loc[, 5] |> p.adjust("BH")

spsf$locEffect <- locEffect
spsf$p.val.adj <- p.val.adj

fill_range <- seq(-3, 8, by = 0.1)

pDetected <- ggplot() +
  geom_sf(data = spsf, aes(color = Nrgn), size = 0.7) +
  scale_color_distiller(palette = "BuPu", direction = 1) +
  geom_point(size = 0.7) +
  labs(color = "Nrgn", title = "Expression of Nrgn\n") +
  theme_void()

pLocEff_direct <- ggplot() +
  geom_sf(data = spsf, aes(color = locEffect), size = 0.7) +
  #geom_sf(data = spsf, aes(color = locEffect), size = 0.3) +
  labs(color = "locI(Nrgn)", title = "Local Moran's I\n(Direct neighbours)") +
  theme_void() +
  scale_color_gradientn(colours = c("blue","white","red"), 
                        values = rescale(c(-1.1,0,4.6)),
                        guide = "colorbar", limits=c(-1.1,4.6))

pLocEff_direct

pLocEff_pVal <- ggplot() +
  geom_sf(data = spsf, aes(color = -log10(p.val.adj)), size = 0.7) +
  scale_color_distiller(palette = "BuGn", direction = 1) +
  labs(color = "-log10(adj. p-val.)", title = "Local Moran's I\n(Adjusted significance levels)") +
  theme_void()

pLocEff_pVal

pLocEff_direct


spsf |>
  ggplot(aes(x = Nrgn, y = locEffect, fill = -log10(p.val.adj))) +
  geom_point() +
  geom_point(colour = "black", pch = 21) +
  scale_fill_distiller(palette = "BuGn", direction = 1) +
  theme_light()



pGeneLoc <- spsf |>
  ggplot(aes(x = Nrgn, y = locEffect, fill = -log10(p.val.adj))) +
  geom_point() +
  geom_point(colour = "black", pch = 21) +
  scale_fill_distiller(palette = "BuGn", direction = 1) +
  theme_light() +
  geom_density_2d(color = "cornflowerblue", bins = 7)


pLocPval <- spsf |>
  ggplot(aes(fill = Nrgn, x = locEffect, y = -log10(p.val.adj))) +
  geom_point() +
  geom_point(colour = "black", pch = 21) +
  scale_fill_distiller(palette = "BuPu", direction = 1) +
  theme_light() +
  geom_density_2d(color = "cornflowerblue", bins = 7)

psupploc <- pGeneLoc + pLocPval


## Plot
p2 <- pDetected + pLocEff_direct + pLocEff_pVal

## Example neighbourhood
## Source: https://pachterlab.github.io/voyager/articles/vig4_cosmx.html
(sfe <- HeNSCLCData())

# Empty cells
colData(sfe)$is_empty <- colData(sfe)$nCounts < 1
# Select, sum negative control probes
(neg_inds <- str_detect(rownames(sfe), "^NegPrb")) %>% sum()


colData(sfe)$prop_neg <- colSums(counts(sfe)[neg_inds, ]) / colData(sfe)$nCounts
# Remove low quality cells
sfe <- sfe[, !sfe$is_empty & sfe$prop_neg < 0.1]
# Re-calculate stats
rowData(sfe)$means <- rowMeans(counts(sfe))
rowData(sfe)$vars <- rowVars(counts(sfe))
rowData(sfe)$is_neg <- neg_inds
# log Counts
sfe <- logNormCounts(sfe)

sfePoly <- sfe
sfeDist <- sfe

colGraph(sfe, "knn10") <- findSpatialNeighbors(sfe,
                                               method = "knearneigh",
                                               dist_type = "idw", k = 10,
                                               style = "W"
)

colGraph(sfeDist, "distance") <- findSpatialNeighbors(sfeDist,
                                                      method = "dnearneigh",
                                                      d1 = 0, d2 = 1000,
                                                      dist_type = "idw",
                                                      style = "W"
)


colGraph(sfePoly, "poly2nb") <- colGeometries(sfePoly)[["cellSeg"]] |>
  poly2nb(snap = 5) |> # snap = 10
  nb2listw(zero.policy = TRUE)


sfe <- runUnivariate(sfe,
                     features = c("KRT17"),
                     type = "localmoran",
                     zero.policy = TRUE
)

bbox_use <- st_bbox(c(xmin = 14000, xmax = 18000, ymin = 158000, ymax = 162000))

# kNN based
(pKnn <- plotLocalResult(sfe, "localmoran",
                         features = c("KRT17"), ncol = 1,
                         colGeometryName = "cellSeg", linewidth = 0.05,
                         bbox = bbox_use
))

pKnn <- pKnn + scale_fill_gradientn(colours = c("#0028A5","#FAFAFA","#FFC845", "#BF0D3E"), 
                                    values = rescale(c(-3,0,5,17)),
                                    guide = "colorbar", limits=c(-3,17)) +
  labs(title = "Local Moran's I\n(10 Nearest Neighbours)", fill = "locI(KRT17)")

pKnn

# poly
sfePoly <- runUnivariate(sfePoly,
                         features = c("KRT17"),
                         type = "localmoran",
                         zero.policy = TRUE
)


(pPoly <- plotLocalResult(sfePoly, "localmoran",
                          features = c("KRT17"), ncol = 1,
                          colGeometryName = "cellSeg", linewidth = 0.05,
                          bbox = bbox_use
))


pPoly <- pPoly + scale_fill_gradientn(colours = c("#0028A5","#FAFAFA","#FFC845", "#BF0D3E"), 
                                      values = rescale(c(-3,0,5,17)),
                                      guide = "colorbar", limits=c(-3,17)) +
  labs(title = "Local Moran's I\n(Contiguos Neighbours)", fill = "locI(KRT17)") 


# distance
sfeDist <- runUnivariate(sfeDist,
                         features = c("KRT17"),
                         type = "localmoran",
                         zero.policy = TRUE
)


(pDist <- plotLocalResult(sfeDist, "localmoran",
                          features = c("KRT17"), ncol = 1,
                          colGeometryName = "cellSeg", linewidth = 0.05,
                          bbox = bbox_use
))


### Figure 2 in main manuscript

pDist <- pDist + scale_fill_gradientn(colours = c("#0028A5","#FAFAFA","#FFC845", "#BF0D3E"), 
                                      values = rescale(c(-3,0,5,17)),
                                      guide = "colorbar", limits=c(-3,17)) +
  labs(title = "Local Moran's I\n(Neighbours in 1000 pixel distance)", fill = "locI(KRT17)") +
  theme_classic()

fig2 <- pDetected + pLocEff_direct + pLocEff_pVal + pPoly + pKnn + pDist +
  plot_layout(nrow = 2, ncol = 3) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag.position = c(0, 1),
    plot.tag = element_text(size = 12, hjust = 0, vjust = 0)
  )

fig2

ggsave("outs/localMoranOverview.pdf", fig2, width = 12, height = 8)


### Supplementary Figure Lattice part
data.frame(
  knn = localResult(sfe, "localmoran", "KRT17")[, "Ii"],
  cont = localResult(sfePoly, "localmoran", "KRT17")[, "Ii"]
) |>
  pivot_longer(everything(), names_to = "neighbourhood", values_to = "Ii") |>
  ggplot(aes(x = Ii)) +
  geom_histogram() +
  facet_grid(neighbourhood ~ .)


logExpHist <- data.frame(KRT17 = logcounts(sfe)["KRT17",]) |> 
  ggplot(aes(x = KRT17)) +
  geom_histogram() +
  theme_light()

# Inspection of distribution of resulting values
localResult(sfe, "localmoran", "KRT17")[, "Ii"] |> 
  round(digits = 5) |> 
  table() |> sort(decreasing = TRUE) |> 
  head(n = 10)

localResult(sfePoly, "localmoran", "KRT17")[, "Ii"] |> 
  round(digits = 5) |> 
  table() |> sort(decreasing = TRUE) |> 
  head(n = 10)

pknncont <- data.frame(
  knn = localResult(sfe, "localmoran", "KRT17")[, "Ii"],
  cont = localResult(sfePoly, "localmoran", "KRT17")[, "Ii"]
) |>
  mutate(zero_neigh = ifelse((cont == 0 & knn != 0), "red", "black")) |> 
  ggplot(aes(x = knn, y = cont, color = zero_neigh)) +
  geom_point(size = 0.5, show.legend = FALSE) +
  geom_abline(intercept = 0, slope = 1, col = "seagreen") +
  theme_light() +
  labs(
    x = "Local Moran's I \n(10 Nearest Neighbours)",
    y = "Local Moran's I \n(Contiguos Neighbours)"
  ) +
  coord_fixed(ratio = 1) +
  geom_density_2d(color = "cornflowerblue") +
  scale_color_manual(values = c("black", "red"))
  
pDistcont <- data.frame(
  dist = localResult(sfeDist, "localmoran", "KRT17")[, "Ii"],
  cont = localResult(sfePoly, "localmoran", "KRT17")[, "Ii"]
) |>
  mutate(zero_neigh = ifelse((cont == 0 & dist != 0), "red", "black")) |> 
  ggplot(aes(x = dist, y = cont, color = zero_neigh)) +
  geom_point(size = 0.5, show.legend = FALSE) +
  geom_abline(intercept = 0, slope = 1, col = "seagreen") +
  geom_density_2d(color = "cornflowerblue") +
  theme_light() +
  labs(
    x = "Local Moran's I \n(Neighbours in 1000 pixel distance)",
    y = "Local Moran's I \n(Contiguos Neighbours)"
  ) +
  coord_fixed(ratio = 1) +
  scale_color_manual(values = c("black", "red"))
  
data.frame(
  knn = localResult(sfe, "localmoran", "KRT17")[, "Ii"],
  cont = localResult(sfePoly, "localmoran", "KRT17")[, "Ii"]
) |>
  ggplot(aes(x = knn, y = cont)) +
  geom_point(size = 0.5, show.legend = FALSE) +
  geom_abline(intercept = 0, slope = 1, col = "cornflowerblue") +
  theme_light() +
  labs(
    x = "Local Moran's I \n(10 Nearest Neighbours)",
    y = "Local Moran's I \n(Contiguos Neighbours)"
  ) +
  coord_fixed(ratio = 1) +
  geom_smooth(method = "lm", col = "red")


# layout for supplementary figure
layout <- 
"
AAAA
AAAA
AAAA
AAAA
AAAA
BBBB
"

# supplementary figure
psupp <- (pGeneLoc + pLocPval + pknncont + pDistcont) /  logExpHist +
  plot_layout(design = layout) +
  plot_annotation(tag_levels = "A")

#plot
psupp

ggsave("outs/localMoranOverviewSupp.pdf", psupp, width = 8, height = 8)

