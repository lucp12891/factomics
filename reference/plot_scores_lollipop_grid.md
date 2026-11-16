# Plot factor scores as lollipop charts for multiple methods

Produces a grid of lollipop plots (one panel per method) showing each
sample's factor score as a segment from a baseline. When `align_to` is
set, each non-reference method displays its **best-matched** factor (by
max \|r\| with the reference) and the sign is aligned automatically.

## Usage

``` r
plot_scores_lollipop_grid(
  bench,
  factor_idx = 1,
  group_col = "Treatment",
  methods = NULL,
  baseline = c("method_mean", "group_mean", "overall_median"),
  palette = NULL,
  legend_position = "right",
  align_to = "auto",
  title_text = "Factor scores",
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

- baseline:

  Baseline: `"method_mean"` (default), `"group_mean"`, or
  `"overall_median"`.

- palette:

  Optional named colour vector for groups.

- legend_position:

  Legend position string (default `"right"`).

- align_to:

  Character. Reference method for factor matching and sign alignment.
  `"auto"` (default) uses MOFA when present, otherwise the first method.
  Set `NULL` to disable.

- title_text:

  Character. Figure title.

- ncol:

  Number of panel columns (default `2`).

## Value

A patchwork ggplot object.
