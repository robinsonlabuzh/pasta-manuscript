# pasta: Pattern Analysis for Spatial Omics daTA

This repository contains code to reproduce analysis and figures of the [pasta preprint](
https://doi.org/10.48550/arXiv.2412.01561).

## Installing necessary R packages
R packages are managed with `renv`. To install R packages you can use `renv::restore()`.

## Data

We use public datasets for all plots. The `STARmap` dataset used in figure 1 can be downloaded from the [singe cell portal](https://singlecell.broadinstitute.org/single_cell/study/SCP1830/spatial-atlas-of-molecular-cell-types-and-aav-accessibility-across-the-whole-mouse-brain) as specified in the [Bansky vignette](https://www.bioconductor.org/packages/release/bioc/vignettes/Banksy/inst/doc/domain-segment.html). All other datasets are downloaded when running `code/utils.R`.

## R Session Info
```
R version 4.5.1 (2025-06-13)
Platform: aarch64-apple-darwin20
Running under: macOS Sequoia 15.5

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: Europe/Zurich
tzcode source: internal

attached base packages:
[1] stats4    stats     graphics  grDevices datasets  utils     methods   base     

other attached packages:
 [1] stringr_1.5.1                   SFEData_1.10.0                  scales_1.4.0                   
 [4] tmap_4.1                        scater_1.36.0                   scuttle_1.18.0                 
 [7] RColorBrewer_1.1-3              viridis_0.6.5                   viridisLite_0.4.2              
[10] ggspavis_1.14.3                 imcdatasets_1.16.0              cytomapper_1.20.0              
[13] EBImage_4.50.0                  sosta_1.0.1                     ggrastr_1.0.2                  
[16] STexampleData_1.16.0            ExperimentHub_2.16.1            AnnotationHub_3.16.1           
[19] BiocFileCache_2.16.1            dbplyr_2.5.0                    patchwork_1.3.1                
[22] SpaNorm_1.2.0                   pals_1.10                       SpatialExperiment_1.18.1       
[25] SingleCellExperiment_1.30.1     SummarizedExperiment_1.38.1     Biobase_2.68.0                 
[28] GenomicRanges_1.60.0            GenomeInfoDb_1.44.1             IRanges_2.42.0                 
[31] S4Vectors_0.46.0                BiocGenerics_0.54.0             generics_0.1.4                 
[34] MatrixGenerics_1.20.0           matrixStats_1.5.0               BiocParallel_1.42.1            
[37] spatialFDA_1.0.0                openxlsx_4.2.8                  spatstat_3.4-0                 
[40] spatstat.linnet_3.3-1           spatstat.model_3.4-0            rpart_4.1.24                   
[43] spatstat.explore_3.5-2          nlme_3.1-168                    spatstat.random_3.4-1          
[46] spatstat.geom_3.5-0             spatstat.univar_3.1-4           spatstat.data_3.1-6            
[49] Voyager_1.9.2                   SpatialFeatureExperiment_1.10.1 ggplot2_3.5.2                  
[52] tidyr_1.3.1                     spdep_1.3-13                    sf_1.0-21                      
[55] spData_2.3.4                    dplyr_1.1.4                         
```
