suppressPackageStartupMessages({library(dplyr)
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
library(patchwork)
library(ggrastr)
})
sfe <- STexampleData::Janesick_breastCancer_Xenium_rep1()

### code adapted from https://lmweber.org/OSTA/pages/crs-spat-stat.html ###
### written by Samuel Gunz and Martin Emons 2025 ###

# load the official 10X annotations
labels <- read.xlsx("https://cdn.10xgenomics.com/raw/upload/v1695234604/Xenium%20Preview%20Data/Cell_Barcode_Type_Matrices.xlsx", sheet = 4)
labels$cell_id <- (labels$Barcode)

# add the cell type labels to the spe
matchedDf <- as.data.frame(colData(sfe)) |>
  left_join(as.data.frame(labels), by = join_by("cell_id"))

colData(sfe) <- DataFrame(matchedDf)

xy <- spatialCoords(sfe)

df <- data.frame(xy, colData(sfe))
df <- df %>% mutate(cell_type = ifelse(Cluster %in% c("DCIS_1", "DCIS_2", "Invasive_Tumor"), Cluster, "Other"))
p0 <- ggplot(df) +
  guides(col=guide_legend(override.aes=list(size=2))) +
  theme(legend.key.size=ggplot2::unit(0, "pt")) +
  labs(col = "cell type", x = "", y = "") +
  theme_light() +
  geom_point(aes(x_centroid, y_centroid, col=cell_type), size=0.3, stroke = 0) +
  coord_fixed() +
  scale_colour_manual(values = c("#0028A5","#FFC845", "#BF0D3E", "grey"),
                      labels = c("DCIS 1", "DCIS 2", "Invasive Tumor", "Other"))


#calculate Lest on the unmarked pattern to compare against
colData(sfe)$unmarked <- "unmarked"
resUnmarked <- calcMetricPerFov(
  sfe,
  selection = "unmarked",
  subsetby = 'sample_id',
  fun = 'Lest',
  marks = 'unmarked',
  rSeq = seq(0, 1000, length.out = 500),
  by = c("sample_id"),
  correction = 'trans'
)

plotMetricPerFov(resUnmarked, theo = TRUE,
                      correction = "trans", x = "r", imageId = "sample_id"
)

##### calculate homo L on the entire tissue ####

resCross <- calcCrossMetricPerFov(
  sfe,
  selection = c("DCIS_1", "DCIS_2", "Invasive_Tumor"),
  subsetby = 'sample_id',
  fun = 'Lcross',
  marks = 'Cluster',
  rSeq = seq(0, 1000, length.out = 500),
  by = c("sample_id"),
  correction = 'trans'
)

resCross$`translation-theoretical` <- resCross$trans - resCross$theo

resCross <- resCross %>% separate(col = selection, into = c("elem1", "elem2"), sep = "and", remove = FALSE)
strip <- ggh4x::strip_themed(background_x = ggh4x::elem_list_rect(fill = c("#0028A5","#FFC845", "#BF0D3E")),
                             text_x = ggh4x::elem_list_text(colour = "white"))


pHomo <- ggplot(resCross, aes(
  x = .data[['r']], y = .data[['translation-theoretical']],
  group = factor(.data[['elem2']])
)) +
  geom_line(aes(color = factor(.data[['elem2']])), linewidth = 1.25) +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(title =
         "Homogeneous L-function in the entire Window",
       colour = "second cell type",
       x = "radius (µm)",
       y = "Translation correction - CSR"
  ) +
  scale_colour_manual(values = c("#0028A5","#FFC845", "#BF0D3E")) +
  # geom_hline(aes(yintercept=0),linetype = "dashed",
  #            color = "black", linewidth = 1) +
  geom_hline(aes(yintercept=0),linetype = "dashed",
             color = "black", linewidth = 1) +
  ggh4x::facet_wrap2(~elem1, strip = strip,) +
  guides(color=guide_legend(override.aes=list(size=2)))

pHomo

pAll <- p0/pHomo
pAll <- pAll + plot_annotation(tag_levels = 'A') +
  plot_layout(heights = c(3, 1))
plot(pAll)

##### calculate inhomo L corrected by the density of the unmarked patter ####

#define intensity of the unmarked pattern for inhomogeneity correction
df <- .speToDf(sfe)
pp <- .dfToppp(df, marks = "Cluster")

#dens_im <- density()
dens_im <- density(unmark(pp))
plot(dens_im)

# set to correct scale as we have less points compared to unmarked
scale <- sum(dens_im$v) / (dens_im$xstep * dens_im$ystep)

resCross <- calcCrossMetricPerFov(
  sfe,
  selection = c("DCIS_1", "DCIS_2", "Invasive_Tumor"),
  subsetby = 'sample_id',
  fun = 'Lcross.inhom',
  marks = 'Cluster',
  rSeq = seq(0, 500, length.out = 500),
  by = c("sample_id"),
  correction = 'trans',
  ncores = 4,
  lambdaX = dens_im * scale
)

resCross$`translation-theoretical` <- resCross$trans - resCross$theo

resCross <- resCross %>% separate(col = selection, into = c("elem1", "elem2"), sep = "and", remove = FALSE)
strip <- ggh4x::strip_themed(background_x = ggh4x::elem_list_rect(fill = c("#0028A5","#FFC845", "#BF0D3E")),
                             text_x = ggh4x::elem_list_text(colour = "white"))


pInhomo <- ggplot(resCross, aes(
  x = .data[['r']], y = .data[['translation-theoretical']],
  group = factor(.data[['elem2']])
)) +
  geom_line(aes(color = factor(.data[['elem2']])), linewidth = 1.25) +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(title =
    "Inhomogeneous L-function in the entire Window",
    colour = "second cell type",
    x = "radius (µm)",
    y = "Translation correction - CSR"
  ) +
  scale_colour_manual(values = c("#0028A5","#FFC845", "#BF0D3E"),
                      labels = c("DCIS 1", "DCIS 2", "Invasive Tumor")) +
  geom_hline(aes(yintercept=0),linetype = "dashed",
             color = "black", linewidth = 1) +
  ggh4x::facet_wrap2(~elem1, strip = strip,) +
  guides(color=guide_legend(override.aes=list(size=2)))

pInhomo

pAll <- p0/pInhomo
pAll <- pAll + plot_annotation(tag_levels = 'A') +
  plot_layout(heights = c(3, 1))
plot(pAll)

##### calculate homo L on the subwindow only ####

## Segment "regions of interest"

library("sosta")

allClusts <- sfe$Cluster |> unique()

imD <- shapeIntensityImage(sfe, marks = "Cluster", imageCol = "sample_id",
                           imageId = "sample01", markSelect = allClusts,
                           bndw = 100)

imD

# take lower threshold than calculated
esThres <- 2.5E-3

# estimated threshold
regions <- reconstructShapeDensitySPE(
  sfe,
  marks = "Cluster",
  imageCol = "sample_id",
  markSelect = allClusts,
  thres = esThres,
  bndw = 100,
  dim = 800
)

plot(regions)

resCross <- calcCrossMetricPerFov(
  sfe,
  selection = c("DCIS_1", "DCIS_2", "Invasive_Tumor"),
  subsetby = 'sample_id',
  fun = 'Lcross',
  marks = 'Cluster',
  rSeq = seq(0, 500, length.out = 500),
  by = c("sample_id"),
  correction = 'trans',
  window = as.owin(regions)
)

resCross$`translation-theoretical` <- resCross$trans - resCross$theo

resCross <- resCross %>% separate(col = selection, into = c("elem1", "elem2"), sep = "and", remove = FALSE)

resCross$elem1  <- resCross$elem1 %>% gsub("_", " ", .)

strip <- ggh4x::strip_themed(background_x = ggh4x::elem_list_rect(fill = c("#0028A5","#FFC845", "#BF0D3E")),
                             text_x = ggh4x::elem_list_text(colour = "white"))

pWindow <- ggplot(resCross, aes(
  x = .data[['r']], y = .data[['translation-theoretical']],
  group = factor(.data[['elem2']])
)) +
  geom_line(aes(color = factor(.data[['elem2']])), linewidth = 1.25) +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(title =
         "Homogeneous L-function in restricted Window",
       colour = "second cell type",
       x = "radius (µm)",
       y = "Translation correction - CSR"
  ) +
  scale_colour_manual(values = c("#0028A5","#FFC845", "#BF0D3E"),
                      labels = c("DCIS 1", "DCIS 2", "Invasive Tumor")) +
  geom_hline(aes(yintercept=0),linetype = "dashed",
             color = "black", linewidth = 1) +
  ggh4x::facet_wrap2(~elem1, strip = strip,) +
  guides(color=guide_legend(override.aes=list(size=2)))

p0 <- p0 + geom_sf(
  data = regions,
  fill = NA,
  color = "black",
  inherit.aes = FALSE, # this is important
  linewidth = 0.8
)

p0 <- rasterize(p0, layers='Point', dpi=200)

pAll <- p0/pWindow
pAll <- pAll + plot_annotation(tag_levels = 'A') +
  plot_layout(heights = c(3, 1))

pAll
ggsave(plot =pAll, "outs/fig3.pdf", width = 10, height = 10, dpi = 100)

pAllOptions <- pHomo / pInhomo / pWindow +
  plot_annotation(tag_levels = 'A')
ggsave(plot =pAllOptions, "outs/fig3complete.pdf", width = 10, height = 15, dpi = 100)
