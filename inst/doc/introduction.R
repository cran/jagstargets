## ---- include = FALSE---------------------------------------------------------
# With the root.dir option below,
# this vignette runs the R code in a temporary directory
# so new files are written to temporary storage
# and not the user's file space.
knitr::opts_knit$set(root.dir = fs::dir_create(tempfile()))
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
skip <- identical(Sys.getenv("NOT_CRAN", unset = "false"), "false") ||
  !requireNamespace("rjags", quietly = TRUE) ||
  !requireNamespace("R2jags", quietly = TRUE)
if (skip) {
  knitr::opts_chunk$set(eval = FALSE)
} else {
  library(R2jags)
}
library(dplyr)
library(targets)
library(jagstargets)

## -----------------------------------------------------------------------------
lines <- "model {
  for (i in 1:n) {
    y[i] ~ dnorm(x[i] * beta, 1)
  }
  beta ~ dnorm(0, 1)
}"
writeLines(lines, "x.jags")

## ---- echo = FALSE------------------------------------------------------------
# Writes the _targets.R file shown in the next code chunk.
library(targets)
tar_script({
  library(jagstargets)
  options(crayon.enabled = FALSE)
  list(
    tar_jags(
      example,
      jags_files = "x.jags",
      parameters.to.save = "beta",
      data = tar_jags_example_data(),
      stdout = R.utils::nullfile(),
      stderr = R.utils::nullfile()
    )
  )
})

## ---- eval = FALSE------------------------------------------------------------
#  # _targets.R
#  library(targets)
#  library(jagstargets)
#  
#  generate_data <- function(n = 10) {
#    true_beta <- stats::rnorm(n = 1, mean = 0, sd = 1)
#    x <- seq(from = -1, to = 1, length.out = n)
#    y <- stats::rnorm(n, x * true_beta, 1)
#    out <- list(n = n, x = x, y = y)
#  }
#  
#  # The _targets.R file ends with a list of target objects
#  # produced by jagstargets::tar_jags(), targets::tar_target(), or similar.
#  list(
#    tar_jags(
#      example,
#      jags_files = "x.jags",
#      parameters.to.save = "beta",
#      data = generate_data()
#    )
#  )

## -----------------------------------------------------------------------------
tar_manifest()

## ---- output = FALSE, message = FALSE-----------------------------------------
tar_visnetwork(targets_only = TRUE)

## ---- output = FALSE----------------------------------------------------------
tar_make()

## -----------------------------------------------------------------------------
tar_read(example_summary_x)

## -----------------------------------------------------------------------------
tar_make()

## -----------------------------------------------------------------------------
write(" ", file = "x.jags", append = TRUE)

## -----------------------------------------------------------------------------
tar_outdated()

## -----------------------------------------------------------------------------
tar_visnetwork(targets_only = TRUE)

## ---- output = FALSE----------------------------------------------------------
tar_make()

## ---- echo = FALSE------------------------------------------------------------
# Writes the _targets.R file shown in the next code chunk.
tar_script({
  library(jagstargets)
  options(crayon.enabled = FALSE)
  tar_option_set(memory = "transient", garbage_collection = TRUE)
  list(
    tar_jags(
      example,
      jags_files = "x.jags",
      parameters.to.save = "beta",
      data = tar_jags_example_data(),
      stdout = R.utils::nullfile(),
      stderr = R.utils::nullfile(),
    ),
    tar_target(
      custom_summary,
      posterior::summarize_draws(
        dplyr::select(example_draws_x, -.draw),
        ~posterior::quantile2(.x, probs = c(0.25, 0.75))
      )
    )
  )
})

## ---- eval = FALSE------------------------------------------------------------
#  # _targets.R
#  library(targets)
#  library(jagstargets)
#  
#  generate_data <- function(n = 10) {
#    true_beta <- stats::rnorm(n = 1, mean = 0, sd = 1)
#    x <- seq(from = -1, to = 1, length.out = n)
#    y <- stats::rnorm(n, x * true_beta, 1)
#    out <- list(n = n, x = x, y = y)
#  }
#  
#  list(
#    tar_jags(
#      example,
#      jags_files = "x.jags",
#      parameters.to.save = "beta",
#      data = generate_data()
#    ),
#    tar_target(
#      custom_summary,
#      posterior::summarize_draws(
#        dplyr::select(example_draws_x, -.draw),
#        ~posterior::quantile2(.x, probs = c(0.25, 0.75))
#      )
#    )
#  )

## -----------------------------------------------------------------------------
tar_visnetwork(targets_only = TRUE)

## ---- output = FALSE, warning = FALSE-----------------------------------------
tar_make()

## -----------------------------------------------------------------------------
tar_read(custom_summary)

## -----------------------------------------------------------------------------
lines <- "model {
  for (i in 1:n) {
    y[i] ~ dnorm(x[i] * x[i] * beta, 1) # Regress on x^2 instead of x.
  }
  beta ~ dnorm(0, 1)
}"
writeLines(lines, "y.jags")

## ---- echo = FALSE------------------------------------------------------------
# Writes the _targets.R file shown in the next code chunk.
tar_script({
  library(targets)
  library(jagstargets)
  list(
    tar_jags(
      example,
      jags_files = c("x.jags", "y.jags"),
      parameters.to.save = "beta",
      data = tar_jags_example_data(),
      stdout = R.utils::nullfile(),
      stderr = R.utils::nullfile(),
    ),
    tar_target(
      custom_summary,
      posterior::summarize_draws(
        dplyr::select(example_draws_x, -.draw),
        ~posterior::quantile2(.x, probs = c(0.25, 0.75))
      )
    )
  )
})

## ---- eval = FALSE------------------------------------------------------------
#  # _targets.R
#  library(targets)
#  library(jagstargets)
#  
#  generate_data <- function(n = 10) {
#    true_beta <- stats::rnorm(n = 1, mean = 0, sd = 1)
#    x <- seq(from = -1, to = 1, length.out = n)
#    y <- stats::rnorm(n, x * true_beta, 1)
#    out <- list(n = n, x = x, y = y)
#  }
#  
#  list(
#    tar_jags(
#      example,
#      jags_files = c("x.jags", "y.jags"),
#      parameters.to.save = "beta",
#      data = generate_data()
#    ),
#    tar_target(
#      custom_summary,
#      posterior::summarize_draws(
#        dplyr::select(example_draws_x, -.draw),
#        ~posterior::quantile2(.x, probs = c(0.25, 0.75))
#      )
#    )
#  )

## -----------------------------------------------------------------------------
tar_visnetwork(targets_only = TRUE)

