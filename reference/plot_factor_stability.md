# Plot the stability of factor methods across seeds

Plot the stability of factor methods across seeds

## Usage

``` r
plot_factor_stability(x, type = c("stability", "recovery"))
```

## Arguments

- x:

  The output of
  [`factor_stability()`](https://lucp12891.github.io/factomics/reference/factor_stability.md).

- type:

  `"stability"` (default): agreement between every pair of seeds.
  `"recovery"`: recovery of the true factors per seed (needs `truth` in
  [`factor_stability()`](https://lucp12891.github.io/factomics/reference/factor_stability.md)).

## Value

A ggplot object.
