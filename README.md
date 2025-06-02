# pasta: Pattern Analysis for Spatial Omics daTA

This repository contains code to reproduce analysis and figures of the [pasta preprint](
https://doi.org/10.48550/arXiv.2412.01561).

## Installing necessary R packages
R packages are managed with `renv`. To install R packages you can use `renv::restore()`.

## R Session Info
```
R version 4.4.1 (2024-06-14)
Platform: aarch64-apple-darwin20
Running under: macOS Sonoma 14.7.1

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: Europe/Zurich
tzcode source: internal

attached base packages:
[1] stats4    stats     graphics  grDevices datasets  utils     methods   base     

other attached packages:
 [1] data.table_1.16.4              Banksy_1.0.0                   imcdatasets_1.12.0            
 [4] cytomapper_1.16.0              sosta_0.99.1                   scales_1.3.0                  
 [7] tidyr_1.3.1                    RColorBrewer_1.1-3             viridis_0.6.5                 
[10] viridisLite_0.4.2              ggspavis_1.10.0                spatialFDA_0.99.1             
[13] MerfishData_1.6.0              EBImage_4.46.0                 magrittr_2.0.3                
[16] stringr_1.5.1                  dixon_0.0-9                    splancs_2.01-45               
[19] sp_2.1-4                       spdep_1.3-8                    spData_2.3.3                  
[22] tmap_3.3-4                     scater_1.32.1                  scran_1.32.0                  
[25] scuttle_1.14.0                 SFEData_1.6.0                  Voyager_1.6.0                 
[28] SpatialFeatureExperiment_1.6.1 ncf_1.3-2                      sf_1.0-19                     
[31] reshape2_1.4.4                 patchwork_1.3.0                STexampleData_1.12.3          
[34] ExperimentHub_2.12.0           AnnotationHub_3.12.0           BiocFileCache_2.12.0          
[37] dbplyr_2.5.0                   RANN_2.6.2                     rlang_1.1.4                   
[40] ggplot2_3.5.1                  dplyr_1.1.4                    mixR_0.2.1                    
[43] spatstat_3.3-0                 spatstat.linnet_3.2-3          spatstat.model_3.3-3          
[46] rpart_4.1.23                   spatstat.explore_3.3-3         nlme_3.1-166                  
[49] spatstat.random_3.3-2          spatstat.geom_3.3-4            spatstat.univar_3.1-1         
[52] spatstat.data_3.1-4            SpatialExperiment_1.14.0       SingleCellExperiment_1.26.0   
[55] SummarizedExperiment_1.34.0    Biobase_2.64.0                 GenomicRanges_1.56.2          
[58] GenomeInfoDb_1.40.1            IRanges_2.38.1                 S4Vectors_0.42.1              
[61] BiocGenerics_0.50.0            MatrixGenerics_1.16.0          matrixStats_1.4.1             
[64] BiocManager_1.30.25           

loaded via a namespace (and not attached):
  [1] R.methodsS3_1.8.2         dichromat_2.0-0.1         tiff_0.1-12              
  [4] urlchecker_1.0.1          goftest_1.2-3             Biostrings_2.72.1        
  [7] HDF5Array_1.32.1          vctrs_0.6.5               digest_0.6.37            
 [10] png_0.1-8                 proxy_0.4-27              pcaPP_2.0-5              
 [13] ggrepel_0.9.6             deldir_2.0-4              renv_1.0.11              
 [16] magick_2.8.5              hdrcde_3.4                MASS_7.3-61              
 [19] httpuv_1.6.15             scico_1.5.0               withr_3.0.2              
 [22] ellipsis_0.3.2            memoise_2.0.1             s2_1.1.7                 
 [25] ggbeeswarm_0.7.2          refund_0.1-37             profvis_0.4.0            
 [28] systemfonts_1.1.0         ragg_1.3.3                lwgeom_0.2-14            
 [31] R.oo_1.27.0               KEGGREST_1.44.1           promises_1.3.2           
 [34] httr_1.4.7                rhdf5filters_1.16.0       ps_1.8.1                 
 [37] rhdf5_2.48.0              rstudioapi_0.17.1         UCSC.utils_1.0.0         
 [40] units_0.8-5               miniUI_0.1.1.1            generics_0.1.3           
 [43] base64enc_0.1-3           processx_3.8.4            stars_0.6-7              
 [46] curl_6.0.1                zlibbioc_1.50.0           ScaledMatrix_1.12.0      
 [49] leafem_0.2.3              polyclip_1.10-7           GenomeInfoDbData_1.2.12  
 [52] SparseArray_1.4.8         fftwtools_0.9-11          xtable_1.8-4             
 [55] desc_1.4.3                pracma_2.4.4              S4Arrays_1.4.1           
 [58] tmaptools_3.1-1           irlba_2.3.5.1             colorspace_2.1-1         
 [61] filelock_1.0.3            isoband_0.2.7             later_1.4.1              
 [64] pbs_1.1                   lattice_0.22-6            mapproj_1.2.11           
 [67] XML_3.99-0.17             svgPanZoom_0.3.4          class_7.3-22             
 [70] pillar_1.10.0             leafsync_0.1.0            compiler_4.4.1           
 [73] beachmat_2.20.0           RSpectra_0.16-2           stringi_1.8.4            
 [76] tensor_1.5                minqa_1.2.8               devtools_2.4.5           
 [79] plyr_1.8.9                fda_6.2.0                 crayon_1.5.3             
 [82] abind_1.4-8               pals_1.9                  locfit_1.5-9.10          
 [85] bit_4.5.0.1               terra_1.8-5               codetools_0.2-20         
 [88] textshaping_0.4.1         BiocSingular_1.20.0       crosstalk_1.2.1          
 [91] leaflet_2.2.2             e1071_1.7-16              fds_1.8                  
 [94] mime_0.12                 splines_4.4.1             Rcpp_1.0.13-1            
 [97] sparseMatrixStats_1.16.0  magic_1.6-1               blob_1.2.4               
[100] BiocVersion_3.19.1        lme4_1.1-35.5             fs_1.6.5                 
[103] nnls_1.6                  DelayedMatrixStats_1.26.0 pkgbuild_1.4.5           
[106] gamm4_0.2-6               RcppHungarian_0.3         tibble_3.2.1             
[109] Matrix_1.7-1              callr_3.7.6               statmod_1.5.0            
[112] svglite_2.1.3             pkgconfig_2.0.3           tools_4.4.1              
[115] aricode_1.0.3             cachem_1.1.0              RSQLite_2.3.9            
[118] DBI_1.2.3                 fastmap_1.2.0             grid_4.4.1               
[121] usethis_3.1.0             shinydashboard_0.7.2      farver_2.1.2             
[124] mgcv_1.9-1                wk_0.9.4                  yaml_2.3.10              
[127] deSolve_1.40              cli_3.6.3                 purrr_1.0.2              
[130] dbscan_1.2-0              lifecycle_1.0.4           uwot_0.2.2               
[133] rainbow_3.8               mvtnorm_1.3-2             bluster_1.14.0           
[136] sessioninfo_1.2.2         DropletUtils_1.24.0       BiocParallel_1.38.0      
[139] gtable_0.3.6              rjson_0.2.23              parallel_4.4.1           
[142] limma_3.60.6              jsonlite_1.8.9            edgeR_4.2.2              
[145] bitops_1.0-9              bit64_4.5.2               spatstat.utils_3.1-1     
[148] BiocNeighbors_1.22.0      ggside_0.3.1              metapod_1.12.0           
[151] dqrng_0.4.1               zeallot_0.1.0             R.utils_2.12.3           
[154] shiny_1.10.0              htmltools_0.5.8.1         rappdirs_0.3.3           
[157] glue_1.8.0                XVector_0.44.0            RCurl_1.98-1.16          
[160] RLRsim_3.1-8              mclust_6.1.1              classInt_0.4-10          
[163] ks_1.14.3                 jpeg_0.1-10               gridExtra_2.3            
[166] sccore_1.0.5              boot_1.3-31               igraph_2.1.2             
[169] R6_2.5.1                  labeling_0.4.3            cluster_2.1.8            
[172] pkgload_1.4.0             Rhdf5lib_1.26.0           memuse_4.2-3             
[175] nloptr_2.1.1              DelayedArray_0.30.1       tidyselect_1.2.1         
[178] vipor_0.4.7               maps_3.4.2.1              raster_3.6-30            
[181] leidenAlg_1.1.4           AnnotationDbi_1.66.0      sfheaders_0.4.4          
[184] rsvd_1.0.5                munsell_0.5.1             KernSmooth_2.23-24       
[187] grpreg_3.5.0              htmlwidgets_1.6.4         spatstat.sparse_3.1-0    
[190] remotes_2.5.0             ggnewscale_0.5.0          beeswarm_0.4.0      
```
