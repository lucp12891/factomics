# Simulate a multi-omics data set with known latent factors

Generates several omics views that share a small number of latent
factors. The first factors separate the sample groups, so a good factor
method should recover the groups, and only a sparse subset of features
in each view carries signal, so feature-level recovery can also be
scored.

## Usage

``` r
simulate_multiomics(
  n_samples = 60,
  views = c(transcriptome = 200, metabolome = 80, microbiome = 50),
  n_groups = 3,
  n_factors = 3,
  n_group_factors = 2,
  effect = 2.5,
  sparsity = 0.2,
  noise = 1,
  seed = 21
)
```

## Arguments

- n_samples:

  Integer. Number of samples (default `60`).

- views:

  Named integer vector giving the number of features per view.

- n_groups:

  Integer. Number of sample groups (default `3`).

- n_factors:

  Integer. Number of true latent factors (default `3`).

- n_group_factors:

  Integer. How many of the factors are shifted by group (default `2`).
  The remaining factors are continuous nuisance variation.

- effect:

  Numeric. Distance of each group centroid from the origin on the group
  factors (default `2.5`).

- sparsity:

  Numeric in (0, 1\]. Fraction of features in each view loaded on each
  factor (default `0.2`).

- noise:

  Numeric. Standard deviation of the Gaussian noise (default `1`).

- seed:

  Integer random seed, or `NULL` (default `21`).

## Value

A list with

- `data`:

  Named list of **features x samples** matrices.

- `metadata`:

  Data frame with `sample`, `group`, `batch` and `age`.

- `truth`:

  List with the true scores `Z`, loadings `W` (per view), and
  `feature_labels`, a named 0/1 vector flagging signal features
  (`"view::feature"` names).

## Examples

``` r
sim <- simulate_multiomics(n_samples = 30, seed = 1)
str(sim$data)
#> List of 3
#>  $ transcriptome: num [1:200, 1:30] -0.155 1.709 1.992 0.316 -3.889 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:200] "gene_001" "gene_002" "gene_003" "gene_004" ...
#>   .. ..$ : chr [1:30] "S001" "S002" "S003" "S004" ...
#>  $ metabolome   : num [1:80, 1:30] -0.513 1.8 -0.924 -1.555 3.274 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:80] "met_001" "met_002" "met_003" "met_004" ...
#>   .. ..$ : chr [1:30] "S001" "S002" "S003" "S004" ...
#>  $ microbiome   : num [1:50, 1:30] -3.149 -0.309 -5.352 -2.73 -0.165 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:50] "otu_001" "otu_002" "otu_003" "otu_004" ...
#>   .. ..$ : chr [1:30] "S001" "S002" "S003" "S004" ...
table(sim$metadata$group)
#> 
#> G1 G2 G3 
#> 10 10 10 
```
