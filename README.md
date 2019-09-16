
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Reproduce Results of Gerard and Hoff (2017)

This folder contains the code necessary to reproduce the results of
Gerard and Hoff (2017). To reproduce these results you will need to:

1.  Download the appropriate R packages.
2.  Run `make`.
3.  Get some coffee.

# Download the appropriate R packages.

You can obtain all of the needed R packages by running the following
code in R:

``` r
install.packages(c("dplyr", "ggplot2", "tidyr", "xtable", 
                   "devtools", "snow", "cate", "ggthemes",
                   "stringr"))
source("https://bioconductor.org/biocLite.R")
biocLite("sva")
devtools::install_github("dcgerard/tensr")
devtools::install_github("dcgerard/hose")
```

To download the code in this repo, click on [this
link](https://github.com/dcgerard/reproduce_sure/archive/master.zip).

# Run `make`

To reproduce all of the results of Gerard and Hoff (2017), simply run
`make` from the terminal (not in the R session). To reproduce the figure
from Section 2, run in the terminal:

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

# Get coffee

Some of the simulations will take awhile to run (2 to 10 hours depending
on how many cores you are using). You should get some coffee\! Here is a
list of some of my favorite places:

  - Chicago
      - [Sawada Coffee](https://www.yelp.com/biz/sawada-coffee-chicago)
      - [Plein Air
        Cafe](https://www.yelp.com/biz/plein-air-cafe-and-eatery-chicago-2)
  - Seattle
      - [Bauhaus
        Ballard](https://www.yelp.com/biz/bauhaus-ballard-seattle)
      - [Cafe Solstice](https://www.yelp.com/biz/cafe-solstice-seattle)
  - Columbus
      - [Yeah, Me Too](https://www.yelp.com/biz/yeah-me-too-columbus)
      - [Stauf’s Coffee
        Roasters](https://www.yelp.com/biz/staufs-coffee-roasters-columbus-2)
      - [Caffe
        Apropos](https://www.yelp.com/biz/caff%C3%A9-apropos-columbus-2)

# Bugs

If you have trouble running this code, then it might be that you need to
update your R packages. I ran these simulations under these settings:

``` r
sessionInfo()
#> R version 3.6.1 (2019-07-05)
#> Platform: x86_64-pc-linux-gnu (64-bit)
#> Running under: Ubuntu 18.04.3 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas/libblas.so.3
#> LAPACK: /usr/lib/x86_64-linux-gnu/libopenblasp-r0.2.20.so
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
#>  [1] stringr_1.4.0       ggthemes_4.2.0      cate_1.1           
#>  [4] sva_3.32.1          BiocParallel_1.18.1 genefilter_1.66.0  
#>  [7] mgcv_1.8-28         nlme_3.1-141        snow_0.4-3         
#> [10] xtable_1.8-4        tidyr_1.0.0         ggplot2_3.2.1      
#> [13] tensr_1.0.1         hose_1.0.0          dplyr_0.8.3        
#> 
#> loaded via a namespace (and not attached):
#>  [1] Rcpp_1.0.2           lattice_0.20-38      esaBcv_1.2.1        
#>  [4] corpcor_1.6.9        assertthat_0.2.1     zeallot_0.1.0       
#>  [7] digest_0.6.20        R6_2.4.0             backports_1.1.4     
#> [10] stats4_3.6.1         RSQLite_2.1.2        evaluate_0.14       
#> [13] pillar_1.4.2         rlang_0.4.0          svd_0.5             
#> [16] lazyeval_0.2.2       annotate_1.62.0      blob_1.2.0          
#> [19] S4Vectors_0.22.1     Matrix_1.2-17        rmarkdown_1.15      
#> [22] splines_3.6.1        RCurl_1.95-4.12      bit_1.1-14          
#> [25] munsell_0.5.0        compiler_3.6.1       xfun_0.9            
#> [28] pkgconfig_2.0.2      BiocGenerics_0.30.0  ruv_0.9.7.1         
#> [31] htmltools_0.3.6      tidyselect_0.2.5     gridExtra_2.3       
#> [34] tibble_2.1.3         IRanges_2.18.2       matrixStats_0.55.0  
#> [37] leapp_1.2            XML_3.98-1.20        crayon_1.3.4        
#> [40] withr_2.1.2          MASS_7.3-51.4        bitops_1.0-6        
#> [43] grid_3.6.1           gtable_0.3.0         lifecycle_0.1.0     
#> [46] DBI_1.0.0            magrittr_1.5         scales_1.0.0        
#> [49] stringi_1.4.3        limma_3.40.6         vctrs_0.2.0         
#> [52] tools_3.6.1          bit64_0.9-7          Biobase_2.44.0      
#> [55] glue_1.3.1           purrr_0.3.2          parallel_3.6.1      
#> [58] survival_2.44-1.1    yaml_2.2.0           AnnotationDbi_1.46.1
#> [61] colorspace_1.4-1     memoise_1.1.0        knitr_1.24
```

As you can see above, I’ve only tried this on Linux.

If you still have difficulty, please submit an
[issue](https://github.com/dcgerard/reproduce_sure/issues).

# References

<div id="refs" class="references">

<div id="ref-gerard2017adaptive">

Gerard, David, and Peter Hoff. 2017. “Adaptive Higher-Order Spectral
Estimators.” *Electron. J. Statist.* 11 (2): 3703–37.
<https://doi.org/10.1214/17-EJS1330>.

</div>

</div>
