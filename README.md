
# <code>Power Analysis</code><br>Power Analysis for Differential Expression in scRNA-seq data

[![](https://img.shields.io/badge/devel%20version-0.0.0.9000-black.svg)](https://github.com/neurogenomics/Power_Analysis_package)
[![](https://img.shields.io/github/languages/code-size/neurogenomics/Power_Analysis_package.svg)](https://github.com/neurogenomics/Power_Analysis_package)
[![](https://img.shields.io/github/last-commit/neurogenomics/Power_Analysis_package.svg)](https://github.com/neurogenomics/Power_Analysis_package/commits/master)
<br> [![R build
status](https://github.com/neurogenomics/Power_Analysis_package/workflows/rworkflows/badge.svg)](https://github.com/neurogenomics/Power_Analysis_package/actions)
<br> [![License: MIT + file
LICENSE](https://img.shields.io/badge/license-MIT%20+%20file%20LICENSE-blue.svg)](https://cran.r-project.org/web/licenses/MIT%20+%20file%20LICENSE)

**Authors:** ***Salman Fawad, Alan Murphy, Yunning Yuan, Hiranyamaya
(Hiru) Dash, Nathan Skene***  
**Updated:** ***Apr-16-2025***

## Introduction

The poweranalysis R package is designed to run robust power analysis for
differential gene expression in scRNA-seq studies and provides tools to
estimate the optimal number of samples and cells needed to achieve
reliable power levels.

Wraps work from
[this](https://github.com/neurogenomics/scRNA-seq_Power_Analysis)
repository into an R package.

To install:

``` r
devtools::install_github("neurogenomics/Power_Analysis")
```

Load:

``` r
library(poweranalysis)
```

## Contact

### [Neurogenomics Lab](https://www.neurogenomics.co.uk)

UK Dementia Research Institute  
Department of Brain Sciences  
Faculty of Medicine  
Imperial College London  
[GitHub](https://github.com/neurogenomics)

<hr>

### Session Info

<details>

``` r
utils::sessionInfo()
```

    ## R version 4.4.3 (2025-02-28)
    ## Platform: aarch64-apple-darwin20
    ## Running under: macOS Sequoia 15.4
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRblas.0.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: Asia/Kolkata
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] gtable_0.3.6        jsonlite_2.0.0      renv_1.1.4         
    ##  [4] dplyr_1.1.4         compiler_4.4.3      BiocManager_1.30.25
    ##  [7] tidyselect_1.2.1    rvcheck_0.2.1       scales_1.3.0       
    ## [10] yaml_2.3.10         fastmap_1.2.0       here_1.0.1         
    ## [13] ggplot2_3.5.1       R6_2.6.1            generics_0.1.3     
    ## [16] knitr_1.50          yulab.utils_0.2.0   tibble_3.2.1       
    ## [19] desc_1.4.3          dlstats_0.1.7       munsell_0.5.1      
    ## [22] rprojroot_2.0.4     pillar_1.10.2       RColorBrewer_1.1-3 
    ## [25] rlang_1.1.5         badger_0.2.4        xfun_0.52          
    ## [28] fs_1.6.5            cli_3.6.4           magrittr_2.0.3     
    ## [31] rworkflows_1.0.6    digest_0.6.37       grid_4.4.3         
    ## [34] rstudioapi_0.17.1   lifecycle_1.0.4     vctrs_0.6.5        
    ## [37] evaluate_1.0.3      glue_1.8.0          data.table_1.17.0  
    ## [40] colorspace_2.1-1    rmarkdown_2.29      tools_4.4.3        
    ## [43] pkgconfig_2.0.3     htmltools_0.5.8.1

</details>
