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

sfe <- sfe[, colSums(counts(sfe)) > 0]
sfe <- sfe[, !is.na(sfe$Cluster)]
sfe <- logNormCounts(sfe)

xy <- spatialCoords(sfe)

df <- data.frame(xy, colData(sfe))
df <- df %>% mutate(cell_type = ifelse(Cluster %in% c("DCIS_1", "DCIS_2", "Invasive_Tumor"), Cluster, "Other"))
p0 <- ggplot(df) +
  guides(col=guide_legend(override.aes=list(size=2))) +
  theme(legend.key.size=ggplot2::unit(0, "pt")) +
  labs(col = "cell type") +
  theme_light() +
  xlab('x') + 
  ylab('y') +
  geom_point(aes(x_centroid, y_centroid, col=cell_type), size=0.3, stroke = 0) +
  coord_fixed() +
  scale_colour_manual(values = c("#0028A5","#FFC845", "#BF0D3E", "grey"))

resCross <- calcCrossMetricPerFov(
  sfe,
  selection = c("DCIS_1", "DCIS_2", "Invasive_Tumor"),
  subsetby = 'sample_id',
  fun = 'Lcross',
  marks = 'Cluster',
  rSeq = seq(0, 1000, length.out = 500),
  by = c("sample_id"),
  correction = 'iso'
)

resCross <- resCross %>% separate(col = selection, into = c("elem1", "elem2"), sep = "and", remove = FALSE)
strip <- ggh4x::strip_themed(background_x = ggh4x::elem_list_rect(fill = c("#0028A5","#FFC845", "#BF0D3E")),
                             text_x = ggh4x::elem_list_text(colour = "white"))

p1 <- ggplot(resCross, aes(
  x = .data[['r']], y = .data[['iso']],
  group = factor(.data[['elem2']])
)) +
  geom_line(aes(color = factor(.data[['elem2']])), linewidth = 1.25) +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(title = paste0(
    resCross$fun, " metric of the pairwise combinations"),
    colour = "second cell type"
  ) +
  scale_colour_manual(values = c("#0028A5","#FFC845", "#BF0D3E")) +
  geom_line(aes(x = .data[['r']], y = theo),
              linetype = "dashed", color = "black", linewidth = 1
  ) + ggh4x::facet_wrap2(~elem1, strip = strip,) +
  guides(color=guide_legend(override.aes=list(size=2))) 

pAll <- p0/p1
pAll <- pAll + plot_annotation(tag_levels = 'A') +
  plot_layout(heights = c(3, 1))

ggsave(plot =pAll, "outs/fig3.pdf", width = 10, height = 10, dpi = 200)

