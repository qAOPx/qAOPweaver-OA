# Direct qKER and ODE replay of the archived demo settings.
# The selected computational functions are kept in upstream.R.

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(2L, 3L)) {
  stop("Usage: Rscript --vanilla r/qker_ode/run.R ARCHIVE_ROOT RESULTS_ROOT [all|qker|ode]",
       call. = FALSE)
}
archive <- normalizePath(args[[1L]], mustWork = TRUE)
results <- args[[2L]]
mode <- if (length(args) == 3L) match.arg(args[[3L]], c("all", "qker", "ode")) else "all"
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
here <- dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE))
source(file.path(here, "upstream.R"), encoding = "UTF-8")

read_input <- function(folder, file) {
  lines <- readLines(file.path(archive, "inputs", folder, file), encoding = "UTF-8")
  lines[[1L]] <- sub(paste0("^", intToUtf8(65279L)), "", lines[[1L]])
  utils::read.csv(text = lines, encoding = "UTF-8", stringsAsFactors = FALSE)
}

case_settings <- jsonlite::read_json(file.path(here, "settings.json"), simplifyVector = TRUE)

write_result <- function(value, ...) {
  path <- file.path(results, ..., "results.json")
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  jsonlite::write_json(value, path, pretty = TRUE, digits = 4)
}

condition_means <- function(rows, expected_rows, observations_per_condition) {
  needed <- c("dataset_id", "ke_id", "dose", "time", "replicate_id", "response")
  if (!all(needed %in% names(rows)) || nrow(rows) != expected_rows ||
      length(unique(rows$dataset_id)) != 1L ||
      any(!is.finite(rows$dose) | !is.finite(rows$time) | !is.finite(rows$response)) ||
      anyDuplicated(rows[c("ke_id", "dose", "time", "replicate_id")])) {
    stop("Expected the complete single-dataset demo input.", call. = FALSE)
  }
  counts <- stats::aggregate(response ~ ke_id + dose + time, rows, length)
  stopifnot(all(counts$response == observations_per_condition))
  means <- stats::aggregate(response ~ ke_id + dose + time, rows, mean)
  means$dataset_id <- rows$dataset_id[[1L]]
  means$time_point <- means$time
  # The source dispatcher's condition identifiers are shared by dose and time.
  # Their literal hashes are immaterial to these joins and are not output.
  means$observation_id <- paste(means$dose, means$time, sep = "@")
  means$replicate_id <- means$observation_id
  means[order(means$ke_id, means$dose, means$time), , drop = FALSE]
}

read_measured <- function() {
  rows <- read_input("measured", "lemna_cccp_ke_response_data_intake_ready.csv")
  condition_means(rows, expected_rows = 360L, observations_per_condition = 5L)
}

read_path <- function(folder, file) {
  path <- read_input(folder, file)
  if (!all(c("ker_id", "from_ke_id", "to_ke_id") %in% names(path)) ||
      nrow(path) != 3L || anyDuplicated(path$ker_id)) {
    stop("Expected the three archived AOP263 edges.", call. = FALSE)
  }
  path$edge_id <- path$ker_id
  path
}

run_qker_case <- function(means) {
  path <- read_path("measured", "lemna_cccp_path_structure_data_intake_ready.csv")
  models <- strsplit(case_settings$qker$models, ",", fixed = TRUE)[[1L]]
  times <- sort(unique(means$time))
  settings <- case_settings$qker$preparation
  reviewed_pairs <- settings$coverage_selection$time_pairs
  stopifnot(
    identical(as.numeric(times), c(2, 4, 7)),
    identical(settings$normalization, "none"),
    identical(settings$replicates, "retain_replicates"),
    setequal(settings$coverage_selection$selected_ker_ids, path$edge_id),
    nrow(reviewed_pairs) == length(times) * nrow(path),
    all(reviewed_pairs$upstream_time == reviewed_pairs$downstream_time)
  )

  combinations <- list()
  for (time in times) {
    fit <- fit_ker(means[means$time == time, , drop = FALSE], path, model_ids = models)
    if (!isTRUE(fit$ok) || !length(fit$combinations)) {
      stop("qKER failed at day ", time, call. = FALSE)
    }
    label <- paste(time, intToUtf8(8594L), time)
    for (combination in fit$combinations) {
      combination$kers <- lapply(combination$kers, function(edge) {
        edge$time_pair <- label
        edge$time_pair_label <- label
        edge
      })
      combination$temporal_pattern <- label
      combination$temporal <- label
      combination$coverage <- 100
      combination$complexity <- combination$kers[[1L]]$model
      combinations[[length(combinations) + 1L]] <- combination
    }
  }
  combinations <- combinations[order(vapply(combinations, function(x) -x$score, numeric(1L)))]
  for (i in seq_along(combinations)) combinations[[i]]$id <- paste0("C", i)
  write_result(list(ok = TRUE, combinations = combinations), "qKER")

  observations <- means
  observations$response_value <- observations$response
  path$ker_label <- path$ker_display_label
  path$from_ke_label <- path$from_ke_name
  path$to_ke_label <- path$to_ke_name
  candidates <- qker_fit_user_data_candidates(observations, path)
  edge_candidates <- lapply(unique(candidates$edge_id), function(edge) {
    rows <- candidates[candidates$edge_id == edge, , drop = FALSE]
    rows <- rows[order(rows$from_time, rows$to_time, rows$candidate_rank), , drop = FALSE]
    list(edge_id = edge, ker_label = rows$ker_label[[1L]], from_ke_id = rows$from_ke_id[[1L]],
         to_ke_id = rows$to_ke_id[[1L]], n_candidates = nrow(rows), candidates = rows)
  })
  write_result(list(ok = TRUE, workflow = "qker_v0320", dataset_id = means$dataset_id[[1L]],
                    candidates = candidates, validation = qker_validate_data_intake_inputs(observations, path),
                    model_candidates_by_edge = edge_candidates),
               "qKER", "all-candidates")
  message("qKER: ", length(combinations), " selected combinations and ", nrow(candidates), " candidate rows.")
}

# The archived single-concentration fits have no forcing schedule. This returns
# the empty-schedule values used by the original derivative evaluator.
ode_sim_forcing_values_at <- function(schedule, time) {
  if (nrow(schedule)) stop("The archived ODE fit has no forcing schedule.", call. = FALSE)
  list(values = numeric(), errors = ode_sim_empty_diag())
}

run_ode_case <- function(kind, means, path) {
  request <- case_settings$ode[[kind]]
  state_ids <- request$states$state_id
  equations <- c(
    KE1446 = "-param_k1 * (KE1446 - param_p1)",
    KE1771 = "-param_k2 * (KE1771 - KE1446)",
    KE1821 = "-param_k3 * (KE1821 - KE1771)",
    KE1521 = "-param_k4 * (KE1521 - KE1821)"
  )
  stopifnot(
    identical(state_ids, names(equations)),
    identical(unlist(request$equations), equations),
    identical(request$graph$edges$edge_id, path$edge_id),
    identical(request$graph$edges$from, path$from_ke_id),
    identical(request$graph$edges$to, path$to_ke_id),
    identical(request$model_family, "linear_qaop"),
    identical(request$mode, "inverse_fitting"),
    identical(request$solver_config$method, "rk4"),
    identical(request$fitting_config$objective, "sum_squared_error"),
    !isTRUE(request$diagnostics_config$profile_likelihood)
  )
  model <- list(equations = data.frame(target_state = state_ids, expression = unname(equations),
                                       is_active = TRUE))
  grid <- request$time_grid
  times <- seq(grid$start, grid$end, by = grid$step)
  initial <- stats::setNames(request$initial_conditions$initial_value, request$initial_conditions$state_id)
  initial <- initial[state_ids]
  parameter_ids <- request$parameters$parameter_id
  bounds <- request$parameter_bounds[match(parameter_ids, request$parameter_bounds$parameter_id), ]
  # The original inverse fitter reads initial_value, then fills missing
  # initials with the midpoint. The request's free-parameter value is separate.
  initial_parameters <- request$parameters$initial_value %||% rep(NA_real_, length(parameter_ids))
  needs_midpoint <- !is.finite(initial_parameters)
  initial_parameters[needs_midpoint] <- ((bounds$lower + bounds$upper) / 2)[needs_midpoint]
  starting <- stats::setNames(initial_parameters, parameter_ids)
  stopifnot(all(request$parameters$parameter_role == "free"), all(is.finite(starting)),
            all(is.finite(initial)), all(is.finite(bounds$lower)), all(is.finite(bounds$upper)),
            all(starting >= bounds$lower & starting <= bounds$upper))

  selected <- means[means$dose == request$dose_value, , drop = FALSE]
  selected <- selected[order(match(selected$ke_id, state_ids), selected$time), , drop = FALSE]
  stopifnot(nrow(selected) == if (kind == "synthetic") 36L else 12L,
            all(selected$time %in% times))
  observed <- data.frame(observation_id = paste0("obs_", seq_len(nrow(selected))),
                         time = selected$time, state_id = selected$ke_id,
                         scenario_id = "default", replicate_id = "",
                         observed_value = selected$response, weight = 1)

  align <- function(parameters) {
    names(parameters) <- names(starting)
    simulation <- ode_sim_integrate(model, parameters, initial, times, method = "rk4")
    if (nrow(simulation$errors)) return(NULL)
    trajectories <- simulation$results
    trajectories$scenario_id <- "default"
    aligned <- ode_fit_align_predictions(trajectories, observed)
    if (nrow(aligned$errors)) return(NULL)
    list(aligned = aligned, trajectories = trajectories)
  }
  objective <- function(parameters) {
    fit <- align(parameters)
    if (is.null(fit) || any(!is.finite(fit$aligned$residuals$residual))) return(1e100)
    sum(fit$aligned$residuals$squared_residual * observed$weight)
  }
  opt <- stats::optim(starting, objective, method = "L-BFGS-B", lower = bounds$lower,
                      upper = bounds$upper, control = list(maxit = request$fitting_config$maxit))
  fit <- align(opt$par)
  if (is.null(fit)) stop("ODE fitting failed for the ", kind, " input.", call. = FALSE)
  eps <- sqrt(.Machine$double.eps)
  fitted_parameters <- data.frame(
    parameter_id = names(starting), estimate = as.numeric(opt$par),
    lower = bounds$lower, upper = bounds$upper, initial_value = as.numeric(starting),
    boundary_hit = abs(opt$par - bounds$lower) <= eps | abs(opt$par - bounds$upper) <= eps
  )
  optimizer <- data.frame(
    optimizer = "stats::optim", method = "L-BFGS-B", convergence_code = opt$convergence,
    convergence_message = opt$message, objective_value = opt$value,
    function_evaluations = unname(opt$counts[["function"]]),
    gradient_evaluations = unname(opt$counts[["gradient"]])
  )
  aligned <- fit$aligned
  fitting_result <- list(
    status = "fitted", fitted_parameters = fitted_parameters,
    optimizer_summary = optimizer, fit_metrics = ode_fit_metrics(aligned$residuals, length(starting), opt$value),
    prediction_alignment = aligned$alignment, fitted_predictions = aligned$predictions,
    residuals = aligned$residuals, trajectories = fit$trajectories
  )
  write_result(list(status = "success", fitting_result = fitting_result), "ODE", kind)
  dir <- file.path(results, "ODE", kind, "Step2")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(data.frame(ke_id = means$ke_id, fit_dose = means$dose,
                              fit_time = means$time, fit_response = means$response),
                   file.path(dir, "preprocessed.csv"), row.names = FALSE)
  message("ODE ", kind, ": ", nrow(observed), " condition means, RMSE ",
          format(fitting_result$fit_metrics$RMSE, digits = 8), ".")
}

measured <- read_measured()
if (mode %in% c("all", "qker")) run_qker_case(measured)
if (mode %in% c("all", "ode")) {
  synthetic <- condition_means(
    read_input("synthetic_ode", "lemna_cccp_ode_step2_demo.csv"),
    expected_rows = 324L, observations_per_condition = 3L
  )
  run_ode_case("synthetic", synthetic,
               read_path("synthetic_ode", "lemna_cccp_ode_path_structure.csv"))
  run_ode_case("measured", measured,
               read_path("measured", "lemna_cccp_path_structure_data_intake_ready.csv"))
}
