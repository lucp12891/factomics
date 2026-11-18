# Settings for ensemble FABIA

Collects the settings of the ensemble FABIA procedure, used by
[`fabia_ensemble()`](https://lucp12891.github.io/factomics/reference/fabia_ensemble.md)
and by `run_factor_pipeline(methods = "FABIA_ENS")`.

## Usage

``` r
fabia_ensemble_control(
  runs = 30,
  sample_fraction = 1,
  similarity = 0.5,
  min_reliability = 0.5,
  membership = 0.5,
  thresZ = 0.5,
  thresL = NULL
)
```

## Arguments

- runs:

  Number of FABIA runs (default `30`). The book used 100.

- sample_fraction:

  Share of samples used in each run. `1` (default) reruns FABIA on all
  samples with a different seed; a value below 1 draws a random
  subsample without replacement for each run.

- similarity:

  Jaccard similarity at which the bicluster tree is cut (default `0.5`:
  biclusters in the same group share at least half of their row-column
  combinations).

- min_reliability:

  Minimum reliability (share of runs in which a group of biclusters was
  found) for a group to be kept (default `0.5`).

- membership:

  Share of a group's biclusters a feature or sample must belong to, to
  be part of the robust bicluster. `1` gives the core (intersection)
  bicluster, a small value such as `0.01` the "super" (union) bicluster,
  and `0.2` the threshold of the biclust package. Default `0.5`.

- thresZ, thresL:

  Thresholds on factor scores and loadings used to turn each FABIA
  factor into a bicluster; passed to
  [`fabia::extractBic()`](https://rdrr.io/pkg/fabia/man/extractBic.html)
  (defaults `0.5` and `NULL`, i.e. FABIA's own choice).

## Value

A named list.

## References

Khamiakova T, Kaiser S, Shkedy Z (2016). Ensemble methods and robust
solutions. In Kasim A, Shkedy Z, Kaiser S, Hochreiter S, Talloen W
(eds), *Applied Biclustering Methods for Big and High-Dimensional Data
Using R*, Chapter 10. Chapman and Hall/CRC.

## Examples

``` r
fabia_ensemble_control(runs = 50, sample_fraction = 0.8)
#> $runs
#> [1] 50
#> 
#> $sample_fraction
#> [1] 0.8
#> 
#> $similarity
#> [1] 0.5
#> 
#> $min_reliability
#> [1] 0.5
#> 
#> $membership
#> [1] 0.5
#> 
#> $thresZ
#> [1] 0.5
#> 
#> $thresL
#> NULL
#> 
```
