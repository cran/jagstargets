# targets::tar_test() runs the test code inside a temporary directory
# to avoid accidentally writing to the user's file space.

# tar_jags() creates a JAGS pipeline that
# runs one MCMC and returns easily consumable outputs,
# such as a data frame of summaries and DIC.
# This test checks that the pipeline is correctly constructed
# and the output is correctly formatted.
targets::tar_test("tar_jags()", {
  skip_on_cran()
  skip_if_not_installed("dplyr")
  skip_if_not_installed("rjags")
  skip_if_not_installed("R2jags")
  tar_jags_example_file(path = "a.jags")
  tar_jags_example_file(path = "b.jags")
  targets::tar_script({
    test_data <- function() {
      out <- tar_jags_example_data(n = 10)
      out$true_beta <- NULL
      out
    }
    list(
      tar_jags(
        model,
        jags_files = c(x = "a.jags", y = "b.jags"),
        parameters.to.save = "beta",
        n.iter = 2e3,
        n.burnin = 1e3,
       data = test_data()
      )
    )
  })
  # Enough targets are created.
  out <- targets::tar_manifest(callr_function = NULL)
  expect_equal(nrow(out), 13L)
  # Nodes in the graph are connected properly.
  out <- targets::tar_network(callr_function = NULL, targets_only = TRUE)$edges
  out <- dplyr::arrange(out, from, to)
  rownames(out) <- NULL
  exp <- tibble::tribble(
    ~from, ~to,
    "model_data", "model_mcmc_x",
    "model_data", "model_mcmc_y",
    "model_file_x", "model_lines_x",
    "model_file_y", "model_lines_y",
    "model_lines_x", "model_mcmc_x",
    "model_lines_y", "model_mcmc_y",
    "model_data", "model_dic_x",
    "model_data", "model_draws_x",
    "model_data", "model_summary_x",
    "model_data", "model_dic_y",
    "model_data", "model_draws_y",
    "model_data", "model_summary_y",
    "model_mcmc_x", "model_dic_x",
    "model_mcmc_x", "model_draws_x",
    "model_mcmc_x", "model_summary_x",
    "model_mcmc_y", "model_dic_y",
    "model_mcmc_y", "model_draws_y",
    "model_mcmc_y", "model_summary_y"
  )
  exp <- dplyr::arrange(exp, from, to)
  rownames(exp) <- NULL
  expect_equal(out, exp)
  # results
  tmp <- capture.output(targets::tar_make(callr_function = NULL))
  expect_equal(targets::tar_read(model_file_x), "a.jags")
  expect_equal(targets::tar_read(model_file_y), "b.jags")
  out <- targets::tar_read(model_data)
  expect_true(is.list(out))
  expect_equal(out$n, 10L)
  expect_equal(length(out$x), 10L)
  expect_equal(length(out$y), 10L)
  expect_true(is.numeric(out$x))
  expect_true(is.numeric(out$y))
  out_x <- targets::tar_read(model_mcmc_x)
  out_y <- targets::tar_read(model_mcmc_y)
  expect_true(inherits(out_x, "rjags"))
  expect_true(inherits(out_y, "rjags"))
  out_x <- targets::tar_read(model_draws_x)
  out_y <- targets::tar_read(model_draws_y)
  expect_true(tibble::is_tibble(out_x))
  expect_true(tibble::is_tibble(out_y))
  expect_equal(nrow(out_x), 3000L)
  expect_equal(nrow(out_y), 3000L)
  expect_true("beta" %in% colnames(out_x))
  expect_true("beta" %in% colnames(out_y))
  out_x <- targets::tar_read(model_summary_x)
  out_y <- targets::tar_read(model_summary_y)
  expect_true(tibble::is_tibble(out_x))
  expect_true(tibble::is_tibble(out_y))
  expect_true("q5" %in% colnames(out_x))
  expect_true("q5" %in% colnames(out_y))
  expect_true(nrow(out_x) < 10)
  expect_true(nrow(out_y) < 10)
  expect_true("mean" %in% colnames(out_x))
  expect_true("mean" %in% colnames(out_y))
  expect_true("beta" %in% out_x$variable)
  expect_true("beta" %in% out_y$variable)
  out_x <- targets::tar_read(model_dic_x)
  out_y <- targets::tar_read(model_dic_y)
  expect_true(tibble::is_tibble(out_x))
  expect_true(tibble::is_tibble(out_y))
  expect_equal(nrow(out_x), 1)
  expect_equal(nrow(out_y), 1)
  expect_true("dic" %in% colnames(out_x))
  expect_true("dic" %in% colnames(out_y))
  # Everything should be up to date.
  expect_equal(targets::tar_outdated(callr_function = NULL), character(0))
  # Change the model. Some targets should invalidate.
  write("", file = "a.jags", append = TRUE)
  out <- targets::tar_outdated(callr_function = NULL)
  exp <- c(
    "model_file_x",
    "model_lines_x",
    "model_dic_x",
    "model_summary_x",
    "model_draws_x",
    "model_mcmc_x"
  )
  expect_equal(sort(out), sort(exp))
  # Change the_data code. Some targets should invalidate.
  targets::tar_script({
    test_data <- function() {
      out <- tar_jags_example_data(n = 10)
      out$true_beta <- NULL
      out
    }
    list(
      tar_jags(
        model,
        jags_files = c(x = "a.jags", y = "b.jags"),
        parameters.to.save = "beta",
        data = test_data(),
        n.chains = 2,
        n.cluster = 2,
        n.iter = 2e3,
        n.burnin = 1e3
      )
    )
  })
  out <- targets::tar_outdated(callr_function = NULL)
  exp <- c(
    exp,
    "model_dic_y",
    "model_summary_y",
    "model_draws_y",
    "model_mcmc_y"
  )
  expect_equal(sort(out), sort(exp))
  # Run with n.cluster > 1.
  capture.output(targets::tar_make(callr_function = NULL))
  expect_true(inherits(targets::tar_read(model_mcmc_y), "rjags"))
})

targets::tar_test("tar_jags() can produced user-defined MCMC summaries", {
  skip_on_cran()
  skip_if_not_installed("dplyr")
  skip_if_not_installed("rjags")
  skip_if_not_installed("R2jags")
  tar_jags_example_file(path = "a.jags")
  tar_jags_example_file(path = "b.jags")
  targets::tar_script({
    test_data <- function() {
      out <- tar_jags_example_data(n = 10)
      out$true_beta <- NULL
      out
    }
    list(
      tar_jags(
        model,
        jags_files = c(x = "a.jags", y = "b.jags"),
        parameters.to.save = "beta",
        n.iter = 2e3,
        n.burnin = 1e3,
        data = test_data(),
        summaries = list(
          custom = ~posterior::quantile2(.x, probs = 0.3),
          custom2 = function(x, my_arg) my_arg
        ),
        summary_args = list(my_arg = 34L)
      )
    )
  })
  capture.output(targets::tar_make(callr_function = NULL))
  out <- targets::tar_read(model_summary_x)
  expect_true("q30" %in% colnames(out))
  expect_true(all(out$custom2 == 34))
})

targets::tar_test("tar_jags() errors is no JAGS file", {
  expect_error(
    tar_jags(
      model,
      jags_files = c(x = "a.jags", y = "b.jags"),
      parameters.to.save = "beta",
      data = test_data(),
      n.chains = 2,
      n.cluster = 2,
      n.iter = 2e3,
      n.burnin = 1e3
    ),
    class = "tar_condition_validate"
  )
})
