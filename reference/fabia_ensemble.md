# Ensemble FABIA: robust biclusters and factors from many FABIA runs

FABIA results depend on the random starting values: two runs on the same
data can return different biclusters. The ensemble approach runs FABIA
many times and keeps only the biclusters found again and again, which
are robust to the starting values (and, with `sample_fraction < 1`, to
the samples used). This is the procedure of Chapter 10 of Kasim et al.
(2016), implemented in the superbiclust and biclust packages, applied to
FABIA.

## Usage

``` r
fabia_ensemble(
  X,
  k = 5,
  control = fabia_ensemble_control(),
  alpha = 0.01,
  cyc = 1000,
  scale_features = TRUE,
  seed = 21,
  verbose = TRUE
)
```

## Arguments

- X:

  Either a single numeric matrix (**features x samples**) or a named
  list of such matrices. All matrices must share the same sample names
  (column names).

- k:

  Number of factors (biclusters) per FABIA run, and the maximum number
  of robust factors returned.

- control:

  Settings, see
  [`fabia_ensemble_control()`](https://lucp12891.github.io/factomics/reference/fabia_ensemble_control.md).

- alpha, cyc:

  FABIA sparseness parameter and number of cycles, as for `"FABIA"` in
  [`run_factor_pipeline()`](https://lucp12891.github.io/factomics/reference/run_factor_pipeline.md).

- scale_features:

  Logical. If `TRUE` (default), each feature is z-scored before fitting.
  Features with zero variance become 0.

- seed:

  Integer. Random seed (default `21`).

- verbose:

  Logical. Print progress messages (default `TRUE`).

## Value

A list with `scores` (samples x factors), `loadings` (features x
factors) and `meta`. `meta` contains:

- `reliability`:

  Reliability of each returned factor (share of runs).

- `biclusters`:

  For each factor, the robust bicluster: its `features`, `samples` and
  `reliability`.

- `groups`:

  All groups of biclusters found: size, number of runs, reliability and
  whether they were kept.

- `similarity`, `tree`:

  The Jaccard similarity matrix of all biclusters and the
  complete-linkage tree (plot it with
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html)).

- `runs`, `control`:

  Settings used.

## Details

Each robust group of biclusters also gives a factor, so the result can
be used like any other factomics method: its loadings and scores are the
sign-aligned averages of the loadings and scores of the FABIA factors in
the group. Factors are ordered by reliability.

This is a shortcut for
`run_factor_pipeline(X, methods = "FABIA_ENS", ...)`.

## References

Khamiakova T, Kaiser S, Shkedy Z (2016). Ensemble methods and robust
solutions. In Kasim A et al. (eds), *Applied Biclustering Methods for
Big and High-Dimensional Data Using R*, Chapter 10. CRC Press.

Hochreiter S et al. (2010). FABIA: factor analysis for bicluster
acquisition. *Bioinformatics* 26(12), 1520-1527.

## See also

[`plot_fabia_ensemble()`](https://lucp12891.github.io/factomics/reference/plot_fabia_ensemble.md),
[`factor_stability()`](https://lucp12891.github.io/factomics/reference/factor_stability.md)

## Examples

``` r
# \donttest{
if (requireNamespace("fabia", quietly = TRUE)) {
  sim <- simulate_multiomics(n_samples = 40, seed = 1)
  fit <- fabia_ensemble(sim$data, k = 3,
                        control = fabia_ensemble_control(runs = 10),
                        cyc = 200, verbose = FALSE)
  fit$meta$reliability
}
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
#>  F1  F2  F3 
#> 1.0 0.9 0.9 
# }
```
