
<!-- README.md is generated from README.Rmd. Please edit that file -->
Reproduce Results of Gerard and Hoff (2015)
===========================================

This folder contains the code necessary to reproduce the results of Gerard and Hoff (2015). To reproduce these results you will need to:

1.  Download the appropriate R packages.
2.  Run `make`.
3.  Get some coffee.

Download the appropriate R packages.
====================================

You can obtain all of the needed R packages by running the following code in R:

``` r
install.packages(c("dplyr", "ggplot2", "tidyr", "xtable", 
                   "devtools", "snow", "cate", "ggthemes",
                   "stringr"))
source("https://bioconductor.org/biocLite.R")
biocLite("sva")
devtools::install_github("dcgerard/tensr")
devtools::install_github("dcgerard/hose")
```

To download the code in this repo, click on [this link](https://github.com/dcgerard/reproduce_sure/archive/master.zip).

Run `make`
==========

To reproduce all of the results of Gerard and Hoff (2015), simply run `make` from the terminal (not in the R session). To reproduce the figure from Section 2, run in the terminal:

``` shell
make change_sv
```

To reproduce the simulation results from the paper, run in the terminal:

``` shell
make sims
```

To reproduce the analysis of NBA statistics, run in the terminal:

``` shell
make nba
```

Get coffee
==========

Some of the simulations will take awhile to run (2 to 10 hours depending on how many cores you are using). You should get some coffee! Here is a list of some of my favorite places:

-   Chicago
    -   [Sawada Coffee](https://www.yelp.com/biz/sawada-coffee-chicago)
    -   [Plein Air Cafe](https://www.yelp.com/biz/plein-air-cafe-and-eatery-chicago-2)
-   Seattle
    -   [Bauhaus Ballard](https://www.yelp.com/biz/bauhaus-ballard-seattle)
    -   [Cafe Solstice](https://www.yelp.com/biz/cafe-solstice-seattle)
-   Columbus
    -   [Yeah, Me Too](https://www.yelp.com/biz/yeah-me-too-columbus)
    -   [Stauf's Coffee Roasters](https://www.yelp.com/biz/staufs-coffee-roasters-columbus-2)
    -   [Caffe Apropos](https://www.yelp.com/biz/caff%C3%A9-apropos-columbus-2)

Bugs
====

If you have trouble running this code, then it might be that you need to update your R packages. I ran these simulations under these settings:

``` r
sessionInfo()
#> R version 3.3.2 (2016-10-31)
#> Platform: x86_64-pc-linux-gnu (64-bit)
#> Running under: Ubuntu 16.04.2 LTS
#> 
#> locale:
#>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
#>  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
#>  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#>  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
#> [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#>  [1] stringr_1.2.0     ggthemes_3.4.0    cate_1.0.4       
#>  [4] sva_3.18.0        genefilter_1.52.1 mgcv_1.8-17      
#>  [7] nlme_3.1-131      snow_0.4-2        xtable_1.8-2     
#> [10] tidyr_0.6.1       ggplot2_2.2.1     tensr_1.0.0      
#> [13] hose_0.1.2        dplyr_0.5.0      
#> 
#> loaded via a namespace (and not attached):
#>  [1] splines_3.3.2        lattice_0.20-34      colorspace_1.3-2    
#>  [4] htmltools_0.3.6      stats4_3.3.2         yaml_2.1.14         
#>  [7] XML_3.98-1.8         survival_2.41-2      rlang_0.1.1         
#> [10] DBI_0.6              BiocGenerics_0.16.1  plyr_1.8.4          
#> [13] munsell_0.4.3        leapp_1.2            gtable_0.2.0        
#> [16] svd_0.4              evaluate_0.10.1      memoise_1.1.0       
#> [19] Biobase_2.30.0       knitr_1.16           IRanges_2.4.8       
#> [22] parallel_3.3.2       AnnotationDbi_1.32.3 esaBcv_1.2.1        
#> [25] Rcpp_0.12.12         corpcor_1.6.8        scales_0.4.1        
#> [28] backports_1.0.5      S4Vectors_0.8.11     annotate_1.48.0     
#> [31] digest_0.6.12        stringi_1.1.2        grid_3.3.2          
#> [34] rprojroot_1.2        tools_3.3.2          magrittr_1.5        
#> [37] lazyeval_0.2.0       tibble_1.3.3         RSQLite_1.1-2       
#> [40] MASS_7.3-45          Matrix_1.2-8         ruv_0.9.6           
#> [43] assertthat_0.2.0     rmarkdown_1.6        R6_2.2.2
```

As you can see above, I've only tried this on Linux.

If you still have difficulty, please submit an [issue](https://github.com/dcgerard/reproduce_sure/issues).

References
==========

Gerard, David, and Peter Hoff. 2015. “Adaptive Higher-Order Spectral Estimators.” *ArXiv Preprint ArXiv:1505.02114*. <http://arxiv.org/abs/1505.02114>.
