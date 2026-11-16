# tests/testthat/helper-data.R
# Shared synthetic multi-omics fixtures.
# Uses <<- / assign(..., envir = .GlobalEnv) so objects are visible in every
# test-file environment regardless of how testthat scopes helper files.

set.seed(42)

.n_samples <<- 12L
.n_feat1   <<- 20L   # microbiome
.n_feat2   <<- 15L   # metabolome
.n_feat3   <<- 10L   # lipidome

.snames <<- paste0("S", seq_len(.n_samples))

.make_mat <- function(nf, tag, snames) {
  matrix(rnorm(nf * length(snames)), nrow = nf,
         dimnames = list(paste0(tag, seq_len(nf)), snames))
}

.data_list <<- list(
  microbiome = .make_mat(.n_feat1, "OTU", .snames),
  metabolome = .make_mat(.n_feat2, "MET", .snames),
  lipidome   = .make_mat(.n_feat3, "LIP", .snames)
)

.single_mat <<- .make_mat(25L, "F", .snames)

.metadata <<- data.frame(
  sample    = .snames,
  Treatment = rep(c("CRC", "CRC+RT", "CRC+IT", "CRC+RT+IT"), each = 3L),
  Batch     = rep(c("A", "B"), times = 6L),
  stringsAsFactors = FALSE
)
