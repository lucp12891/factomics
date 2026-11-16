# List the factor-analysis methods supported by factomics

List the factor-analysis methods supported by factomics

## Usage

``` r
list_factor_methods()
```

## Value

A data frame with one row per method: its name (as passed to `methods`
in
[`run_factor_pipeline()`](https://lucp12891.github.io/factomics/reference/run_factor_pipeline.md)),
the R package that implements it, its family, a short description and
whether that package is installed.

## Examples

``` r
list_factor_methods()
#> +----------------------------+                                          
#> |............................|                                          
#> |............................|                                          
#> |..............########......|  #######    #    ######    ###      #    
#> |..............########......|  #         # #   #     #    #      # #   
#> |.....####.....########......|  #        #   #  #     #    #     #   #  
#> |.....####.....########......|  #####   #     # ######     #    #     # 
#> |.....####...................|  #       ####### #     #    #    ####### 
#> |.....####...........###.....|  #       #     # #     #    #    #     # 
#> |....................###.....|  #       #     # ######    ###   #     # 
#> |....................###.....|                                          
#> |............................|                                          
#> +----------------------------+                                          
#> Citation: S. Hochreiter et al.,
#> FABIA: Factor Analysis for Bicluster Acquisition,
#> Bioinformatics 26(12):1520-1527, 2010.
#> BibTex: enter 'toBibtex(citation("fabia"))'
#> 
#> Homepage: http://www.bioinf.jku.at/software/fabia/fabia.html
#> 
#> FABIA Package Version 2.58.0
#>   method    package                  family
#> 1    PCA      stats    Matrix factorisation
#> 2    MFA FactoMineR    Matrix factorisation
#> 3  FABIA      fabia Bayesian / biclustering
#> 4    GFA        GFA                Bayesian
#> 5   MOFA      MOFA2                Bayesian
#> 6  SGCCA   mixOmics       Correlation-based
#> 7   JIVE     r.jive    Matrix factorisation
#>                                                description installed
#> 1                 Block-weighted PCA on concatenated views      TRUE
#> 2              Multiple Factor Analysis (Escofier & Pages)      TRUE
#> 3                Factor Analysis for Bicluster Acquisition      TRUE
#> 4 Group Factor Analysis (Bayesian, view-specific sparsity)      TRUE
#> 5    Multi-Omics Factor Analysis v2 (needs Python mofapy2)      TRUE
#> 6        Sparse Generalised Canonical Correlation Analysis      TRUE
#> 7                 Joint and Individual Variation Explained      TRUE
```
