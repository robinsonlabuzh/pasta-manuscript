# pasta: Pattern Analysis for Spatial Omics daTA

This repository contains code to reproduce analysis and figures of the [pasta preprint](
https://doi.org/10.48550/arXiv.2412.01561).

## Installing necessary R packages
R packages are managed with `renv`. To install R packages you can use `renv::restore()`.

## Data

We use public datasets for all plots. The `STARmap` dataset used in figure 1 can be downloaded from the [singe cell portal](https://singlecell.broadinstitute.org/single_cell/study/SCP1830/spatial-atlas-of-molecular-cell-types-and-aav-accessibility-across-the-whole-mouse-brain) as specified in the [Bansky vignette](https://www.bioconductor.org/packages/release/bioc/vignettes/Banksy/inst/doc/domain-segment.html). All other datasets are downloaded when running `code/utils.R`.

## R Session Info
```
R Under development (unstable) (2025-01-30 r87669)
Platform: aarch64-apple-darwin20
Running under: macOS Sequoia 15.5

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: Europe/Zurich
tzcode source: internal

attached base packages:
[1] stats4    stats     graphics  grDevices datasets  utils     methods   base     

other attached packages:
 [1] SpaNorm_1.2.0                  pals_1.10                      sosta_1.0.0                   
 [4] STexampleData_1.12.3           ExperimentHub_2.16.0           AnnotationHub_3.16.0          
 [7] BiocFileCache_2.16.0           dbplyr_2.5.0                   ggrastr_1.0.2                 
[10] patchwork_1.3.0                SpatialExperiment_1.18.1       SingleCellExperiment_1.30.1   
[13] SummarizedExperiment_1.38.1    Biobase_2.68.0                 GenomicRanges_1.60.0          
[16] GenomeInfoDb_1.44.0            IRanges_2.42.0                 S4Vectors_0.46.0              
[19] BiocGenerics_0.54.0            generics_0.1.4                 MatrixGenerics_1.20.0         
[22] matrixStats_1.5.0              BiocParallel_1.42.0            spatialFDA_1.0.0              
[25] openxlsx_4.2.8                 spatstat_3.3-3                 spatstat.linnet_3.2-6         
[28] spatstat.model_3.3-6           rpart_4.1.24                   spatstat.explore_3.4-3        
[31] nlme_3.1-167                   spatstat.random_3.4-1          spatstat.geom_3.4-1           
[34] spatstat.univar_3.1-3          spatstat.data_3.1-6            Voyager_1.10.0                
[37] SpatialFeatureExperiment_1.9.7 ggplot2_3.5.2                  tidyr_1.3.1                   
[40] spdep_1.3-11                   sf_1.0-21                      spData_2.3.4                  
[43] dplyr_1.1.4                   
```
