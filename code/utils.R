suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(spatstat)
  library(mixR)
  library(dplyr)
  library(ggplot2)
  library(rlang)
  #library(seg)
  library(stats)
  library(RANN)
  library(STexampleData)
  library(patchwork)
  library(reshape2)
  library(sf)
  library(ncf)
  #library(rgeoda)
  library(Voyager)
  library(SpatialFeatureExperiment)
  library(SFEData)
  library(scran)
  library(scater)
  library(tmap)
  library(spdep)
  library(dixon)
  library(stringr)
  library(magrittr)
})

# create imaging SPE
eh <- ExperimentHub()
q <- query(eh, "MERFISH")
df <- eh[["EH7546"]]

i <- seq_len(9)
cd <- data.frame(df[, i], row.names = 1)

# set sample identifiers
id <- grep("Bregma", names(cd))
names(cd)[id] <- "sample_id"

# rename spatial coordinates
xy <- grep("Centroid", names(cd))
xy <- names(cd)[xy] <- c("x", "y")

# simplify annotations
cd$cluster_id <- cd$Cell_class
for (. in c("Endothelial", "OD Mature", "OD Immature"))
  cd$cluster_id[grep(., cd$cluster_id)] <- .

# extract & sparsify assay data
y <- data.frame(df[, -i], row.names = df[, 1])
y <- as(t(as.matrix(y)), "dgCMatrix")

# construct SPE
(spe <- SpatialExperiment(
  assays = list(exprs  = y),
  spatialCoordsNames = xy,
  colData = cd))

# Define the directory and file paths
data_path <- "data"
out_path <- "outs"

# Check if the directory exists, and create it if it doesn't
if (!dir.exists(data_path)) {
  dir.create(data_path)
}

if (!dir.exists(out_path)) {
  dir.create(out_path)
}

# Save RDS
saveRDS(spe, "data/spe.rds")

# Define helper functions
.ppp <- \(spe, marks = NULL) {
  xy <- spatialCoords(spe)
  w <- owin(range(xy[, 1]), range(xy[, 2]))
  m <- if (is.character(marks)) {
    stopifnot(
      length(marks) == 1,
      marks %in% c(
        gs <- rownames(spe),
        cd <- names(colData(spe))))
    if (marks %in% gs) {
      assay(spe)[marks, ]
    } else if (marks %in% cd) {
      spe[[marks]]
    } else stop("'marks' should be in ",
      "'rownames(.)' or 'names(colData(.))'")
  }
  ppp(xy[, 1], xy[, 2], window = w, marks = factor(m))
}


pValuesHotspotMarks <- function(pp, alpha = 0.05){
  # Code source: https://idblr.rbind.io/post/pvalues-spatial-segregation/

  # Significant p-values assumming normality of the Poisson process
  ## relrisk() computes standard errors based on asymptotic theory, assuming a Poisson process
  # call relative risk function
  f1 <- relrisk(pp_sel,se=TRUE)

  z <- qnorm(alpha/2, lower.tail = F)     # z-statistic
  f1$u <- f1$estimate + z*f1$SE           # Upper CIs
  f1$l <- f1$estimate - z*f1$SE           # Lower CIs
  mu_0 <- as.vector(table(spatstat.geom::marks(pp))/pp$n) # null expectations by type
  f1$p <- f1$estimate # copy structure of pixels, replace values
  for (i in 1:length(f1$p)) {
    f1$p[[i]]$v <- factor(ifelse(mu_0[i] > f1$u[[i]]$v, "lower",
                                ifelse( mu_0[i] < f1$l[[i]]$v, "higher", "none")),
                          levels = c("lower", "none", "higher"))
  }
  return(f1)
}

pValuesHotspot <- function(pp, alpha = 0.05){
  # Code source: https://idblr.rbind.io/post/pvalues-spatial-segregation/
  # density estimate for all marks
  f1 <- density(unmark(pp), se = TRUE)
  # Significant p-values assumming normality of the Poisson process
  z <- qnorm(alpha/2, lower.tail = F)     # z-statistic
  f1$u <- f1$estimate + z*f1$SE           # Upper CIs
  f1$l <- f1$estimate - z*f1$SE           # Lower CIs
  mu_0 <- intensity(unmark(pp)) # null expectations by type
  f1$p <- f1$estimate # copy structure of pixels, replace values
  f1$p$v <- factor(ifelse(mu_0 > f1$u$v, "lower",
                          ifelse( mu_0 < f1$l$v, "higher", "none")),
                   levels = c("lower", "none", "higher"))

  return(f1)
}
