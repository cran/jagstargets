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
    y[i] ~ dnorm(beta[1] + x[i] * beta[2], 1)
  }
  for (i in 1:2) {
    beta[i] ~ dnorm(0, 1)
  }
}"
writeLines(lines, "model.jags")

## ---- echo = FALSE------------------------------------------------------------
# Writes the _targets.R file shown in the next code chunk.
library(targets)
tar_script({
  library(jagstargets)
  options(crayon.enabled = FALSE)
  tar_option_set(memory = "transient", garbage_collection = TRUE)
  generate_data <- function(n = 10L) {
    beta <- stats::rnorm(n = 2, mean = 0, sd = 1)
    x <- seq(from = -1, to = 1, length.out = n)
    y <- stats::rnorm(n, beta[1] + x * beta[2], 1)
    .join_data <- list(beta = beta)
    list(n = n, x = x, y = y, .join_data = .join_data)
  }
  list(
    tar_jags_rep_summary(
      model,
      "model.jags",
      data = generate_data(),
      parameters.to.save = "beta",
      batches = 5, # Number of branch targets.
      reps = 4, # Number of model reps per branch target.
      stdout = R.utils::nullfile(),
      stderr = R.utils::nullfile(),
      variables = "beta",
      summaries = list(
        ~posterior::quantile2(.x, probs = c(0.025, 0.975))
      )
    )
  )
})

## ---- eval = FALSE------------------------------------------------------------
#  # _targets.R
#  library(targets)
#  library(jagstargets)
#  options(crayon.enabled = FALSE)
#  # Use computer memory more sparingly:
#  tar_option_set(memory = "transient", garbage_collection = TRUE)
#  
#  generate_data <- function(n = 10L) {
#    beta <- stats::rnorm(n = 2, mean = 0, sd = 1)
#    x <- seq(from = -1, to = 1, length.out = n)
#    y <- stats::rnorm(n, beta[1] + x * beta[2], 1)
#    # Elements of .join_data get joined on to the .join_data column
#    # in the summary output next to the model parameters
#    # with the same names.
#    .join_data <- list(beta = beta)
#    list(n = n, x = x, y = y, .join_data = .join_data)
#  }
#  
#  list(
#    tar_jags_rep_summary(
#      model,
#      "model.jags",
#      data = generate_data(),
#      parameters.to.save = "beta",
#      batches = 5, # Number of branch targets.
#      reps = 4, # Number of model reps per branch target.
#      variables = "beta",
#      summaries = list(
#        ~posterior::quantile2(.x, probs = c(0.025, 0.975))
#      )
#    )
#  )

## -----------------------------------------------------------------------------
tar_visnetwork()

## ---- output = FALSE, warning = FALSE-----------------------------------------
tar_make()

## -----------------------------------------------------------------------------
tar_load(model)
model

## -----------------------------------------------------------------------------
library(dplyr)
model %>%
  group_by(variable) %>%
  dplyr::summarize(coverage = mean(q2.5 < .join_data & .join_data < q97.5))

## ---- echo = FALSE------------------------------------------------------------
# Writes the _targets.R file shown in the next code chunk.
library(targets)
tar_script({
  library(jagstargets)
  options(crayon.enabled = FALSE)
  tar_option_set(
    packages = "dplyr",
    memory = "transient",
    garbage_collection = TRUE
  )
  generate_data <- function(n = 10L) {
    beta <- stats::rnorm(n = 2, mean = 0, sd = 1)
    x <- seq(from = -1, to = 1, length.out = n)
    y <- stats::rnorm(n, beta[1] + x * beta[2], 1)
    # Elements of .join_data get joined on to the .join_data column
    # in the summary output next to the model parameters
    # with the same names.
    .join_data <- list(beta = beta)
    list(n = n, x = x, y = y, .join_data = .join_data)
  }
  list(
    tar_jags_rep_summary(
      model,
      "model.jags",
      data = generate_data(),
      parameters.to.save = "beta",
      batches = 5, # Number of branch targets.
      reps = 4, # Number of model reps per branch target.
      stdout = R.utils::nullfile(),
      stderr = R.utils::nullfile(),
      variables = "beta",
      summaries = list(
        ~posterior::quantile2(.x, probs = c(0.025, 0.975))
      )
    ),
    tar_target(
      coverage,
      model %>%
        group_by(variable) %>%
        summarize(
          coverage = mean(q2.5 < .join_data & .join_data < q97.5),
          .groups = "drop"
        )
    )
  )
})

## ---- eval = FALSE------------------------------------------------------------
#  # _targets.R
#  # packages needed to define the pipeline:
#  library(targets)
#  library(jagstargets)
#  
#  tar_option_set(
#    packages = "dplyr", # packages needed to run the pipeline
#    memory = "transient", # memory efficiency
#    garbage_collection = TRUE # memory efficiency
#  )
#  
#  generate_data <- function(n = 10L) {
#    beta <- stats::rnorm(n = 2, mean = 0, sd = 1)
#    x <- seq(from = -1, to = 1, length.out = n)
#    y <- stats::rnorm(n, beta[1] + x * beta[2], 1)
#    # Elements of .join_data get joined on to the .join_data column
#    # in the summary output next to the model parameters
#    # with the same names.
#    .join_data <- list(beta = beta)
#    list(n = n, x = x, y = y, .join_data = .join_data)
#  }
#  
#  list(
#    tar_jags_rep_summary(
#      model,
#      "model.jags",
#      data = generate_data(),
#      parameters.to.save = "beta",
#      batches = 5, # Number of branch targets.
#      reps = 4, # Number of model reps per branch target.
#      variables = "beta",
#      summaries = list(
#        ~posterior::quantile2(.x, probs = c(0.025, 0.975))
#      )
#    ),
#    tar_target(
#      coverage,
#      model %>%
#        group_by(variable) %>%
#        summarize(
#          coverage = mean(q2.5 < .join_data & .join_data < q97.5),
#          .groups = "drop"
#        )
#    )
#  )

## -----------------------------------------------------------------------------
tar_visnetwork()

## ---- output = FALSE, warning = FALSE-----------------------------------------
tar_make()

## -----------------------------------------------------------------------------
tar_read(coverage)

## -----------------------------------------------------------------------------
lines <- "model {
  for (i in 1:n) {
    y[i] ~ dnorm(beta[1] + x[i] * x[i] * beta[2], 1) # Regress on x^2, not x.
  }
  for (i in 1:2) {
    beta[i] ~ dnorm(0, 1)
  }
}"
writeLines(lines, "model2.jags")

## ---- echo = FALSE------------------------------------------------------------
# Writes the _targets.R file shown in the next code chunk.
library(targets)
tar_script({
  library(jagstargets)
  options(crayon.enabled = FALSE)
  tar_option_set(
    packages = "dplyr",
    memory = "transient",
    garbage_collection = TRUE
  )
  generate_data <- function(n = 10L) {
    beta <- stats::rnorm(n = 2, mean = 0, sd = 1)
    x <- seq(from = -1, to = 1, length.out = n)
    y <- stats::rnorm(n, beta[1] + x * beta[2], 1)
    # Elements of .join_data get joined on to the .join_data column
    # in the summary output next to the model parameters
    # with the same names.
    .join_data <- list(beta = beta)
    list(n = n, x = x, y = y, .join_data = .join_data)
  }
  list(
    tar_jags_rep_summary(
      model,
      c("model.jags", "model2.jags"), # another model
      data = generate_data(),
      parameters.to.save = "beta",
      batches = 5,
      reps = 4,
      stdout = R.utils::nullfile(),
      stderr = R.utils::nullfile(),
      variables = "beta",
      summaries = list(
        ~posterior::quantile2(.x, probs = c(0.025, 0.975))
      )
    ),
    tar_target(
      coverage,
      model %>%
        group_by(.name) %>%
        summarize(coverage = mean(q2.5 < .join_data & .join_data < q97.5))
    )
  )
})

## ---- eval = FALSE------------------------------------------------------------
#  # _targets.R
#  # packages needed to define the pipeline:
#  library(targets)
#  library(jagstargets)
#  
#  tar_option_set(
#    packages = "dplyr", # packages needed to run the pipeline
#    memory = "transient", # memory efficiency
#    garbage_collection = TRUE # memory efficiency
#  )
#  
#  generate_data <- function(n = 10L) {
#    beta <- stats::rnorm(n = 2, mean = 0, sd = 1)
#    x <- seq(from = -1, to = 1, length.out = n)
#    y <- stats::rnorm(n, beta[1] + x * beta[2], 1)
#    # Elements of .join_data get joined on to the .join_data column
#    # in the summary output next to the model parameters
#    # with the same names.
#    .join_data <- list(beta = beta)
#    list(n = n, x = x, y = y, .join_data = .join_data)
#  }
#  
#  list(
#    tar_jags_rep_summary(
#      model,
#      c("model.jags", "model2.jags"), # another model
#      data = generate_data(),
#      parameters.to.save = "beta",
#      batches = 5,
#      reps = 4,
#      variables = "beta",
#      summaries = list(
#        ~posterior::quantile2(.x, probs = c(0.025, 0.975))
#      )
#    ),
#    tar_target(
#      coverage,
#      model %>%
#        group_by(.name) %>%
#        summarize(coverage = mean(q2.5 < .join_data & .join_data < q97.5))
#    )
#  )

## -----------------------------------------------------------------------------
tar_visnetwork()

