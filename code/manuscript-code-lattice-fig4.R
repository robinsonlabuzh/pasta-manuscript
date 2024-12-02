library(spatialFDA)
library(sosta)
library(ggplot2)
library(spatstat)
library(dplyr)
library(patchwork)
library(sf)


# source: http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
cbPalette <- c("#56B4E9", "#F0E442","#CC79A7", "#E69F00", "#009E73", "#0072B2", "#D55E00", "#999999")

spe <- imcdatasets::Damond_2019_Pancreas("spe", full_dataset = FALSE)
spe

islet <- reconstructShapeDensityImage(
  spe,
  marks = "cell_category",
  image_col = "image_name",
  image_id = "E25",
  mark_select = "islet"
) 

spe_sel <- spe[, spe[["image_name"]] %in% c("E25")]
spe_sel_df <- cbind(spatialCoords(spe_sel), colData(spe_sel)) |> as.data.frame()


spe_sel_df$cell_category %>% unique()
# color palette
group.colors <- c(exocrine = "#56B4E9", immune = "#E69F00",
                  other ="#F0E442", unknown = "#999999", islet = "#CC79A7")

p_eq <- ggplot() +
  geom_point(data = spe_sel_df,
             aes(x = cell_x, y = cell_y, color = cell_category), size = 0.75) +
  theme_light() +
  scale_colour_manual(values=group.colors) +
  theme_void() + labs(color = "") +
  theme(legend.position="bottom") +
  coord_equal()


st_fov <- st_bbox(c(xmin = min(spe_sel_df$cell_x),
                    xmax = max(spe_sel_df$cell_x), 
                    ymin = min(spe_sel_df$cell_y),
                    ymax = max(spe_sel_df$cell_y))) |> st_as_sfc()

p_cell_cat <- ggplot() +
  geom_point(data = spe_sel_df,
             aes(x = cell_x, y = cell_y, color = cell_category), size = 0.75) +
  theme_light() +
  scale_colour_manual(values=group.colors) +
  geom_sf(data = st_fov, fill = NA, color = "#D55E00", linewidth = 1, linetype = "dashed") +
  theme_classic() + labs(color = "", x = "", y = "") +
  theme(legend.position="bottom")

p_cell_cat

idx <- st_intersects(st_as_sf(spe_sel_df, coords = c("cell_x", "cell_y")),islet, sparse = FALSE)


p_cell_islet <- ggplot() +
  geom_point(data = spe_sel_df[idx,],
             aes(x = cell_x, y = cell_y, color = cell_category), size = 0.75) +
  theme_light() +
  scale_colour_manual(values=group.colors) +
  geom_sf(data = islet, fill = NA, color = "#0072B2", linewidth = 1, linetype = "dashed") +
  theme_classic() + labs(color = "", x = "", y = "") +
  theme(legend.position="bottom") +
  scale_x_continuous(limits = c(min(spe_sel_df$cell_x), max(spe_sel_df$cell_x))) +
  scale_y_continuous(limits = c(min(spe_sel_df$cell_y), max(spe_sel_df$cell_y)))



p_cell_islet

p_cell_cat + p_cell_islet  +   plot_annotation(tag_levels = c('A')) +
  #plot_layout(guides = "collect") &
  theme(legend.position='bottom') +
  theme(legend.text=element_text(size=12.5))

##

p_eq | p_cell_cat

# Option 1: Homogeneous global
spe_sub <- subset(spe, , image_name == "E25")

df_sub <- .speToDf(spe_sub)
pp <- .dfToppp(df_sub, marks = "cell_category")

pp_sub <- subset(pp, marks == "islet")

plot(pp)
plot(pp_sub)

plot(Kest(pp_sub))

metric_res1 <- calcMetricPerFov(spe_sub, c("islet"),
                               subsetby = "image_number",
                               fun = "Kest", marks = "cell_category",
                               rSeq = seq(0, 120, length.out = 100), 
                               by = c("patient_stage", "patient_id", "image_number"),
                               ncores = 1
)

p1.1 <- plotMetricPerFov(metric_res1,
                         correction = "iso", x = "r",
                         imageId = "image_number",
                         theo = TRUE) + scale_color_manual(values = "#D55E00") 

p1.1[["layers"]][[2]]$aes_params$linewidth <- 1
p1.1[["layers"]][[1]]$aes_params$linewidth <- 1

p1.1

p1.2 <- plotMetricPerFov(metric_res1,
                         correction = "border", x = "r",
                         imageId = "image_number",
                         theo = TRUE)  + scale_color_manual(values = "#D55E00")

p1.2[["layers"]][[2]]$aes_params$linewidth <- 1
p1.2[["layers"]][[1]]$aes_params$linewidth <- 1


p1.1 + p1.2

p <- ggplot(iris, aes(x = Sepal.Length, y = Sepal.Width)) +
  geom_line(color = "blue")

p

p + geom_line(linewidth = 2)


# Option 2: Inhomogeneous global
# density on all cells
plot(Kinhom(pp_sub, lambda = density(pp, bw.diggle(pp))))

# density on subset
plot(density(pp_sub, bw.diggle(pp_sub)))
plot(Kinhom(pp_sub, lambda = density(pp_sub, bw.diggle(pp_sub))))
plot(Kinhom(unmark(pp_sub), correction = "best"))


metric_res2 <- calcMetricPerFov(spe_sub, "islet",
                                subsetby = "image_number",
                                fun = "Kinhom", marks = "cell_category",
                                rSeq = seq(0, 120, length.out = 100),
                                by = c("patient_stage", "patient_id", "image_number"),
                                ncores = 1
)

### 

p2.1 <- plotMetricPerFov(metric_res2, correction = "iso", 
                          x = "r", imageId = "image_number",
                          theo = TRUE) + scale_color_manual(values = "#D55E00")
p2.2 <- plotMetricPerFov(metric_res2, correction = "trans", 
                       x = "r", imageId = "image_number",
                       theo = TRUE) + scale_color_manual(values = "#D55E00")
p2.3 <- plotMetricPerFov(metric_res2, correction = "border", 
                         x = "r", imageId = "image_number",
                         theo = TRUE) + scale_color_manual(values = "#D55E00")
p2.4 <- plotMetricPerFov(metric_res2, correction = "bord.modif", 
                         x = "r", imageId = "image_number",
                         theo = TRUE) + scale_color_manual(values = "#D55E00")


## 

p2.1[["layers"]][[2]]$aes_params$linewidth <- 1
p2.1[["layers"]][[1]]$aes_params$linewidth <- 1

p2.2[["layers"]][[2]]$aes_params$linewidth <- 1
p2.2[["layers"]][[1]]$aes_params$linewidth <- 1

p2.3[["layers"]][[2]]$aes_params$linewidth <- 1
p2.3[["layers"]][[1]]$aes_params$linewidth <- 1

p2.4[["layers"]][[2]]$aes_params$linewidth <- 1
p2.4[["layers"]][[1]]$aes_params$linewidth <- 1


(p2.1+p2.2)/(p2.3+p2.4)


# Option 3: Set window
plot(islet)

spe_sub <- subset(spe, , image_name == "E25")
df_sub <- .speToDf(spe_sub)

pp <- .dfToppp(df_sub, marks = "cell_category")
pp_sub <- subset(pp, marks == "islet", drop = TRUE)

# Set custom window and exclude all cells outside
pp_win <- as.ppp(superimpose(pp_sub, W =  as.owin(islet)))
plot(pp_win)
plot(Kest(pp_win))
plot(Kinhom(pp_win))


metric_res3 <- .extractMetric(.speToDf(spe_sub), c("islet"),
                             fun = "Kest", marks = "cell_category",
                             rSeq = seq(0, 70, length.out = 100),
                             by = c("patient_stage", "patient_id", "image_number"),
                             window = as.owin(islet))

metric_res3 <- calcMetricPerFov(spe_sub, "islet",
                                subsetby = "image_number",
                                fun = "Kest", marks = "cell_category",
                                rSeq = seq(0, 70, length.out = 100),
                                by = c("patient_stage", "patient_id", "image_number"),
                                ncores = 1,
                                window = as.owin(islet)
)

metric_res4 <- calcMetricPerFov(spe_sub, "islet",
                                subsetby = "image_number",
                                fun = "Kinhom", marks = "cell_category",
                                rSeq = seq(0, 70, length.out = 100),
                                by = c("patient_stage", "patient_id", "image_number"),
                                ncores = 1,
                                window = as.owin(islet)
)



p3.1 <- plotMetricPerFov(metric_res3, correction = "iso", 
                 x = "r", imageId = "image_number",
                 theo = TRUE) + scale_color_manual(values = "#0072B2")

p4.1 <- plotMetricPerFov(metric_res4, correction = "iso", 
                         x = "r", imageId = "image_number",
                         theo = TRUE) + scale_color_manual(values = "#0072B2")
p4.2 <- plotMetricPerFov(metric_res4, correction = "trans", 
                         x = "r", imageId = "image_number",
                         theo = TRUE) + scale_color_manual(values = "#0072B2")
p4.3 <- plotMetricPerFov(metric_res4, correction = "border", 
                         x = "r", imageId = "image_number",
                         theo = TRUE) + scale_color_manual(values = "#0072B2")
p4.4 <- plotMetricPerFov(metric_res4, correction = "bord.modif", 
                         x = "r", imageId = "image_number",
                         theo = TRUE) + scale_color_manual(values = "#0072B2")

p3.1[["layers"]][[2]]$aes_params$linewidth <- 1
p3.1[["layers"]][[1]]$aes_params$linewidth <- 1

p4.1[["layers"]][[2]]$aes_params$linewidth <- 1
p4.1[["layers"]][[1]]$aes_params$linewidth <- 1

p4.2[["layers"]][[2]]$aes_params$linewidth <- 1
p4.2[["layers"]][[1]]$aes_params$linewidth <- 1

p4.3[["layers"]][[2]]$aes_params$linewidth <- 1
p4.3[["layers"]][[1]]$aes_params$linewidth <- 1

p4.4[["layers"]][[2]]$aes_params$linewidth <- 1
p4.4[["layers"]][[1]]$aes_params$linewidth <- 1


p3.1 + p4.1

(p4.1 + p4.2) / (p4.3 + p4.4)

# All of them in one
p4_top <- p_cell_cat + (p1.1 + p1.2) / (p2.1 + p2.2)
p4_bot <- p_cell_islet + (p3.1 + p4.1)
p4_tot_v1 <- p4_top / p4_bot
p4_tot_v1


layout <- 
"
ABC
"

p4_up <- (p_cell_cat  + p1.1 + p2.1) + plot_layout(design = layout)
p4_down <- (p_cell_islet + p3.1 + p4.1) + plot_layout(design = layout)


# Plot for paper
p4_tot <- p4_up / p4_down
p4_tot <- p4_tot + plot_annotation(tag_levels = "A")

p4_tot

ggsave("outs/fig4.pdf", p4_tot, width = 10, height = 7)


# Supplementary
p4_supp <- (p_cell_cat + (p2.1+p2.2)/(p2.3+p2.4)) /
  (p_cell_islet + (p4.1 + p4.2) / (p4.3 + p4.4)) + plot_annotation(tag_levels = "A")

p4_supp
ggsave("outs/fig4-supp.pdf", p4_supp, width = 12, height = 10)
  
