# Adjusted Rand index between two partitions

Adjusted Rand index between two partitions

## Usage

``` r
adjusted_rand_index(a, b)
```

## Arguments

- a, b:

  Vectors of cluster labels of equal length.

## Value

A single number; 1 for identical partitions, about 0 for random.

## Examples

``` r
adjusted_rand_index(c(1, 1, 2, 2), c("a", "a", "b", "b"))
#> [1] 1
```
