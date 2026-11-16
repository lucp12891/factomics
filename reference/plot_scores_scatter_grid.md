# Plot factor scores for multiple methods (scatter grid)

Produces a grid of scatter plots (one panel per method) showing factor
scores per sample, coloured by a metadata grouping variable. When
`align_to` is set, each non-reference method automatically displays its
**matched** factor — i.e. the factor most correlated (by \|r\|) with
factor `factor_idx` of the reference method — and its sign is flipped if
the raw correlation is negative. Use `align_to = NULL` to disable
matching and plot each method's factor `factor_idx` as-is.

## Usage

``` r
plot_scores_scatter_grid(
  bench,
  factor_idx = 1,
  group_col = "Treatment",
  methods = NULL,
  panel_tags = NULL,
  tag_position = c("tl_in", "none"),
  title_text = "Factor scores",
  palette = NULL,
  legend_position = "right",
  align_to = "auto",
  ncol = NULL
)
```

## Arguments

- bench:

  Named list returned by
  [`attach_scores_metadata()`](https://lucp12891.github.io/factomics/reference/attach_scores_metadata.md).

- factor_idx:

  Integer. Factor index to display for the **reference** method
  (`align_to`). Other methods will show their best-matching factor.

- group_col:

  Character. Column from `scores_with_meta` for colouring.

- methods:

  Character vector of method names to plot. `NULL` (default) uses every
  method that has `scores_with_meta`.

- panel_tags:

  Optional character vector of panel labels (one per method).

- tag_position:

  `"tl_in"` (default) or `"none"`.

- title_text:

  Character. Figure title.

- palette:

  Optional named colour vector for groups.

- legend_position:

  Legend position string (default `"right"`).

- align_to:

  Character. Reference method for factor matching and sign alignment.
  `"auto"` (default) uses MOFA when present, otherwise the first method.
  Set `NULL` to disable.

- ncol:

  Number of panel columns (default `2`).

## Value

A patchwork ggplot object.
