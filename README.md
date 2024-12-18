# pasta: Pattern Analysis for Spatial Omics daTA

This repository contains code to reproduce analysis and figures for the [pasta preprint](
https://doi.org/10.48550/arXiv.2412.01561).

## Installing necessary R packages

- `spatialFDA` has to be installed via devtools::install_github("mjemons/spatialFDA").
- `sosta` has to be installed via devtools::install_github("sgunz/sosta").

## R Session Info
```
R version 4.3.1 (2023-06-16)
Platform: aarch64-apple-darwin20 (64-bit)
Running under: macOS Sonoma 14.7.1

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: Europe/Zurich
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices datasets  utils     methods  
[7] base     

loaded via a namespace (and not attached):
  [1] DBI_1.1.3                     bitops_1.0-7                 
  [3] gridExtra_2.3                 rlang_1.1.1                  
  [5] magrittr_2.0.3                scater_1.28.0                
  [7] matrixStats_1.0.0             compiler_4.3.1               
  [9] RSQLite_2.3.1                 DelayedMatrixStats_1.22.6    
 [11] png_0.1-8                     vctrs_0.6.4                  
 [13] reshape2_1.4.4                stringr_1.5.0                
 [15] pkgconfig_2.0.3               crayon_1.5.2                 
 [17] fastmap_1.1.1                 dbplyr_2.3.4                 
 [19] XVector_0.40.0                ellipsis_0.3.2               
 [21] scuttle_1.10.3                utf8_1.2.3                   
 [23] promises_1.2.1                ggbeeswarm_0.7.2             
 [25] bit_4.0.5                     zlibbioc_1.46.0              
 [27] cachem_1.0.8                  beachmat_2.16.0              
 [29] GenomeInfoDb_1.36.4           blob_1.2.4                   
 [31] later_1.3.1                   DelayedArray_0.26.7          
 [33] BiocParallel_1.34.2           interactiveDisplayBase_1.38.0
 [35] irlba_2.3.5.1                 parallel_4.3.1               
 [37] SFEData_1.2.0                 R6_2.5.1                     
 [39] stringi_1.7.12                GenomicRanges_1.52.1         
 [41] Rcpp_1.0.11                   SummarizedExperiment_1.30.2  
 [43] IRanges_2.34.1                httpuv_1.6.11                
 [45] Matrix_1.5-4.1                tidyselect_1.2.0             
 [47] rstudioapi_0.15.0             abind_1.4-5                  
 [49] yaml_2.3.7                    viridis_0.6.4                
 [51] codetools_0.2-19              curl_5.1.0                   
 [53] lattice_0.21-8                tibble_3.2.1                 
 [55] plyr_1.8.9                    Biobase_2.60.0               
 [57] shiny_1.7.5.1                 KEGGREST_1.40.1              
 [59] BiocFileCache_2.8.0           ExperimentHub_2.8.1          
 [61] Biostrings_2.68.1             pillar_1.9.0                 
 [63] BiocManager_1.30.22           filelock_1.0.2               
 [65] MatrixGenerics_1.12.3         renv_1.0.3                   
 [67] stats4_4.3.1                  generics_0.1.3               
 [69] RCurl_1.98-1.12               BiocVersion_3.17.1           
 [71] S4Vectors_0.38.2              ggplot2_3.4.4                
 [73] sparseMatrixStats_1.12.2      munsell_0.5.0                
 [75] scales_1.3.0                  xtable_1.8-4                 
 [77] glue_1.6.2                    tools_4.3.1                  
 [79] AnnotationHub_3.8.0           BiocNeighbors_1.18.0         
 [81] ScaledMatrix_1.8.1            grid_4.3.1                   
 [83] AnnotationDbi_1.62.2          colorspace_2.1-0             
 [85] SingleCellExperiment_1.22.0   GenomeInfoDbData_1.2.10      
 [87] beeswarm_0.4.0                BiocSingular_1.16.0          
 [89] vipor_0.4.5                   cli_3.6.1                    
 [91] rsvd_1.0.5                    rappdirs_0.3.3               
 [93] fansi_1.0.5                   S4Arrays_1.0.6               
 [95] viridisLite_0.4.2             dplyr_1.1.3                  
 [97] gtable_0.3.4                  digest_0.6.33                
 [99] BiocGenerics_0.46.0           ggrepel_0.9.4                
[101] memoise_2.0.1                 htmltools_0.5.6.1            
[103] lifecycle_1.0.3               httr_1.4.7                   
[105] mime_0.12                     bit64_4.0.5
```
