# Default control settings for RABFA

Default control settings for RABFA

## Usage

``` r
rabfa_control(
  nu1 = 0,
  nu2 = 1,
  eta = 10,
  eps = 0.2,
  sigma_m2 = 100,
  phi = 1,
  mh_sd = 1,
  thin = 1
)
```

## Arguments

- nu1:

  Prior mean of the log shrinkage parameters; `exp(nu1)` is the typical
  Laplace rate. Larger values shrink loadings more (default `0`).

- nu2:

  Prior variance scale of the log shrinkage parameters. It sets how much
  each loading's shrinkage can adapt to the data (default `1`).

- eta:

  Concentration of the graph prior: confidence in the graph (default
  `10`, as in the paper).

- eps:

  Controls the strength of the prior correlation between neighbouring
  shrinkage parameters (default `0.2`, as in the paper).

- sigma_m2:

  Prior variance of the feature means (default `100`).

- phi:

  Prior parameter for the noise precisions (default `1`).

- mh_sd:

  Standard deviation of the random-walk proposal for the log shrinkage
  parameters (default `1`).

- thin:

  Keep every `thin`-th draw after burn-in (default `1`).

## Value

A named list.
