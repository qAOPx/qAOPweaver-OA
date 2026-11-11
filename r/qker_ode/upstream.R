# Selected function definitions from qAOPweaver.
# Source revision 0143dd56dd21a5268b17ecbe8383647b1c6323c9

# qAOPweaver/qAOPweaverR/R/report_contract_engine_v0368.R lines 1-3
`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || (length(x) == 1L && is.na(x))) y else x
}

# qAOPweaver/qAOPweaverR/R/frontend_contract_v0380.R lines 8-8
qaop_v0380_num <- function(x) suppressWarnings(as.numeric(x))

# qAOPweaver/qAOPweaverR/R/frontend_contract_v0380.R lines 10-14
qaop_v0380_finite_or_null <- function(x, digits = 6) {
  x <- suppressWarnings(as.numeric(x))[1]
  if (!is.finite(x)) return(NULL)
  round(x, digits)
}

# qAOPweaver/qAOPweaverR/R/frontend_contract_v0380.R lines 24-27
qaop_v0380_first_col <- function(x, choices) {
  hit <- intersect(choices, names(x))
  if (length(hit)) hit[1] else NULL
}

# qAOPweaver/qAOPweaverR/R/frontend_contract_v0380.R lines 29-53
qaop_v0380_normalize_ke_data <- function(ke_data) {
  if (!is.data.frame(ke_data)) stop("ke_data must be a data.frame.", call. = FALSE)
  out <- ke_data
  aliases <- list(
    ke_id = c("ke_id", "endpoint_id", "node_id", "ke_internal_id"),
    dose = c("dose", "exposure_value", "dose_value", "concentration"),
    response = c("response", "response_value", "value", "analysis_response"),
    time = c("time", "time_point", "timepoint", "timepoint_day"),
    replicate = c("replicate", "replicate_id", "rep")
  )
  for (nm in names(aliases)) {
    if (!nm %in% names(out)) {
      col <- qaop_v0380_first_col(out, aliases[[nm]])
      if (!is.null(col)) out[[nm]] <- out[[col]]
    }
  }
  if (!"dataset_id" %in% names(out)) out$dataset_id <- "frontend_dataset"
  if (!"observation_id" %in% names(out)) out$observation_id <- seq_len(nrow(out))
  if (!"time" %in% names(out)) out$time <- 0
  if (!"replicate" %in% names(out)) out$replicate <- seq_len(nrow(out))
  out$dose <- qaop_v0380_num(out$dose)
  out$response <- qaop_v0380_num(out$response)
  out$time <- qaop_v0380_num(out$time)
  out
}

# qAOPweaver/qAOPweaverR/R/frontend_contract_v0380.R lines 55-71
qaop_v0380_normalize_path_data <- function(path_data) {
  if (!is.data.frame(path_data)) stop("path_data must be a data.frame.", call. = FALSE)
  out <- path_data
  aliases <- list(
    from_ke_id = c("from_ke_id", "from", "from_node_id", "upstream_ke_id"),
    to_ke_id = c("to_ke_id", "to", "to_node_id", "downstream_ke_id")
  )
  for (nm in names(aliases)) {
    if (!nm %in% names(out)) {
      col <- qaop_v0380_first_col(out, aliases[[nm]])
      if (!is.null(col)) out[[nm]] <- out[[col]]
    }
  }
  if (!"edge_id" %in% names(out)) out$edge_id <- paste0("E", seq_len(nrow(out)))
  if (!"dataset_id" %in% names(out)) out$dataset_id <- "frontend_dataset"
  out
}

# qAOPweaver/qAOPweaverR/R/frontend_contract_v0380.R lines 180-188
qaop_v0380_model_map <- function(model_id) {
  switch(model_id,
    "LL.4" = "linear_lm",
    "W1.4" = "log_linear_lm",
    "BC.5" = "quadratic_lm",
    "LN.4" = "log_linear_lm",
    model_id
  )
}

# qAOPweaver/qAOPweaverR/R/frontend_contract_v0380.R lines 190-214
qaop_v0380_fit_one <- function(dat, model_id) {
  dat <- dat[stats::complete.cases(dat[, c("x", "y")]), , drop = FALSE]
  if (nrow(dat) < 3) stop("Not enough paired rows for model fitting.", call. = FALSE)
  internal <- qaop_v0380_model_map(model_id)
  fit <- if (identical(internal, "linear_lm")) {
    stats::lm(y ~ x, data = dat)
  } else if (identical(internal, "log_linear_lm")) {
    stats::lm(y ~ log1p(x), data = dat)
  } else if (identical(internal, "quadratic_lm")) {
    stats::lm(y ~ x + I(x^2), data = dat)
  } else {
    stop("Unsupported qKER model_id: ", model_id, call. = FALSE)
  }
  s <- suppressWarnings(summary(fit))
  cf <- stats::coef(fit)
  eq <- if (identical(internal, "linear_lm")) {
    sprintf("y = %.4f + %.4f x", cf[1], cf[2])
  } else if (identical(internal, "log_linear_lm")) {
    sprintf("y = %.4f + %.4f log1p(x)", cf[1], cf[2])
  } else {
    sprintf("y = %.4f + %.4f x + %.4f x^2", cf[1], cf[2], cf[3])
  }
  pval <- if (!is.null(s$coefficients) && nrow(s$coefficients) >= 2 && ncol(s$coefficients) >= 4) s$coefficients[2, 4] else NA_real_
  list(eq = eq, r2 = s$r.squared %||% NA_real_, aic = stats::AIC(fit), pval = pval, n = nrow(dat))
}

# qAOPweaver/qAOPweaverR/R/frontend_contract_v0380.R lines 216-234
qaop_v0380_pair_edge <- function(data, edge, ke_col, dose_col, response_col) {
  up <- data[as.character(data[[ke_col]]) == as.character(edge$from_ke_id), , drop = FALSE]
  down <- data[as.character(data[[ke_col]]) == as.character(edge$to_ke_id), , drop = FALSE]
  keys <- intersect(c("dataset_id", dose_col, "time", "replicate"), intersect(names(up), names(down)))
  if (!length(keys)) keys <- dose_col
  xu <- up[, unique(c(keys, response_col)), drop = FALSE]
  yd <- down[, unique(c(keys, response_col)), drop = FALSE]
  names(xu)[names(xu) == response_col] <- "x"
  names(yd)[names(yd) == response_col] <- "y"
  pair <- merge(xu, yd, by = keys, sort = FALSE)
  if (!nrow(pair) && dose_col %in% names(up) && dose_col %in% names(down)) {
    xu <- stats::aggregate(up[response_col], up[dose_col], mean, na.rm = TRUE)
    yd <- stats::aggregate(down[response_col], down[dose_col], mean, na.rm = TRUE)
    names(xu)[names(xu) == response_col] <- "x"
    names(yd)[names(yd) == response_col] <- "y"
    pair <- merge(xu, yd, by = dose_col, sort = FALSE)
  }
  pair
}

# qAOPweaver/qAOPweaverR/R/frontend_contract_v0380.R lines 239-305
fit_ker <- function(
  data,
  path,
  model_ids = c("LL.4", "W1.4", "BC.5", "LN.4"),
  ke_col = "ke_id",
  dose_col = "dose",
  response_col = "response",
  ...
) {
  tryCatch({
    dat <- qaop_v0380_normalize_ke_data(data)
    path <- qaop_v0380_normalize_path_data(path)
    combos <- list()
    for (mi in seq_along(model_ids)) {
      mid <- model_ids[[mi]]
      kers <- list()
      notes <- character()
      for (ei in seq_len(nrow(path))) {
        edge <- path[ei, , drop = FALSE]
        pair <- qaop_v0380_pair_edge(dat, edge, ke_col, dose_col, response_col)
        fitted <- tryCatch(qaop_v0380_fit_one(pair, mid), error = function(e) e)
        if (inherits(fitted, "error")) {
          notes <- c(notes, paste(edge$edge_id, conditionMessage(fitted), sep = ": "))
          next
        }
        kers[[length(kers) + 1L]] <- list(
          id = paste0("KER-", length(kers) + 1L),
          edge = as.character(edge$edge_id),
          from = as.character(edge$from_ke_id),
          to = as.character(edge$to_ke_id),
          model = as.character(mid),
          eq = fitted$eq,
          r2 = qaop_v0380_finite_or_null(fitted$r2),
          aic = qaop_v0380_finite_or_null(fitted$aic),
          pval = if (is.finite(fitted$pval)) {
            if (fitted$pval < 0.001) "<0.001" else format(round(fitted$pval, 4), nsmall = 4)
          } else NULL,
          n = as.integer(fitted$n)
        )
      }
      if (length(kers)) {
        r2s <- vapply(kers, function(k) as.numeric(k$r2 %||% NA_real_), numeric(1))
        aics <- vapply(kers, function(k) as.numeric(k$aic %||% NA_real_), numeric(1))
        mean_r2 <- mean(r2s, na.rm = TRUE)
        mean_aic <- mean(aics, na.rm = TRUE)
        score <- mean_r2 - (rank(mean_aic, ties.method = "first") * 0.001)
        combos[[length(combos) + 1L]] <- list(
          id = paste0("C", length(combos) + 1L),
          score = qaop_v0380_finite_or_null(score, 3),
          r2 = qaop_v0380_finite_or_null(mean_r2, 3),
          aic = qaop_v0380_finite_or_null(mean_aic, 3),
          unc = if (is.finite(mean_r2) && mean_r2 >= 0.8) "Low" else if (is.finite(mean_r2) && mean_r2 >= 0.5) "Medium" else "High",
          note = if (length(notes)) paste(c("Fit completed with skipped edge/model attempts.", notes), collapse = " ") else "Fit completed for available edges.",
          kers = kers
        )
      }
    }
    if (length(combos)) {
      ord <- order(vapply(combos, function(x) as.numeric(x$score %||% -Inf), numeric(1)), decreasing = TRUE)
      combos <- combos[ord]
      for (i in seq_along(combos)) combos[[i]]$id <- paste0("C", i)
    }
    list(ok = TRUE, combinations = combos)
  }, error = function(e) {
    list(ok = FALSE, combinations = list(list(id = "C1", score = NULL, r2 = NULL, aic = NULL, unc = "High", note = conditionMessage(e), kers = list())))
  })
}

# qAOPweaver/qAOPweaverR/R/qker_user_data_workflow.R lines 1-52
qker_validate_data_intake_inputs <- function(ke_response_data, path_structure_data) {
  ke_required <- c("dataset_id", "observation_id", "dose", "ke_id", "response_value")
  path_required <- c("dataset_id", "edge_id", "from_ke_id", "to_ke_id")
  ke_names <- names(ke_response_data)
  path_names <- names(path_structure_data)
  ke_ok <- all(ke_required %in% ke_names)
  path_ok <- all(path_required %in% path_names)
  dataset_id_consistent <- ke_ok && path_ok &&
    length(unique(ke_response_data$dataset_id)) == 1 &&
    length(unique(path_structure_data$dataset_id)) == 1 &&
    identical(unique(as.character(ke_response_data$dataset_id)), unique(as.character(path_structure_data$dataset_id)))
  path_nodes <- if (path_ok && "ke_id" %in% ke_names) unique(c(path_structure_data$from_ke_id, path_structure_data$to_ke_id)) else character()
  path_nodes_present <- length(path_nodes) > 0 && all(path_nodes %in% unique(ke_response_data$ke_id))
  response_numeric <- "response_value" %in% ke_names && is.numeric(ke_response_data$response_value)
  dose_valid <- "dose" %in% ke_names && is.numeric(ke_response_data$dose) && all(is.finite(ke_response_data$dose))
  time_valid <- !("time_point" %in% ke_names) || (is.numeric(ke_response_data$time_point) && all(is.finite(ke_response_data$time_point)))
  unsupported_structure <- FALSE
  if (path_ok) {
    nodes <- unique(c(as.character(path_structure_data$from_ke_id), as.character(path_structure_data$to_ke_id)))
    children <- split(as.character(path_structure_data$to_ke_id), as.character(path_structure_data$from_ke_id))
    visited <- setNames(rep(FALSE, length(nodes)), nodes)
    active <- visited
    visit <- function(node) {
      if (isTRUE(active[[node]])) return(TRUE)
      if (isTRUE(visited[[node]])) return(FALSE)
      visited[[node]] <<- TRUE
      active[[node]] <<- TRUE
      for (child in children[[node]]) {
        if (!is.null(child) && visit(child)) return(TRUE)
      }
      active[[node]] <<- FALSE
      FALSE
    }
    unsupported_structure <- any(vapply(nodes, visit, logical(1)))
  }
  validation_pass <- ke_ok && path_ok && dataset_id_consistent && path_nodes_present && response_numeric && dose_valid && time_valid && !unsupported_structure
  list(
    required_ke_columns_present = ke_ok,
    required_path_columns_present = path_ok,
    dataset_id_consistent = dataset_id_consistent,
    path_nodes_present_in_ke_data = path_nodes_present,
    response_value_numeric = response_numeric,
    dose_valid = dose_valid,
    observation_id_present = "observation_id" %in% ke_names,
    time_point_handling_valid = time_valid,
    missing_values_reported = anyNA(ke_response_data) || anyNA(path_structure_data),
    duplicates_reported = any(duplicated(ke_response_data[, intersect(c("dataset_id", "observation_id", "ke_id", "time_point", "replicate_id"), ke_names), drop = FALSE])),
    unsupported_structure_reported_without_crash = unsupported_structure,
    validation_pass = validation_pass,
    notes = if (validation_pass) "Data Intake-compatible qKER inputs validated." else "Input validation failed; see boolean columns."
  )
}

# qAOPweaver/qAOPweaverR/R/qker_user_data_workflow.R lines 54-56
qker_v0320_mgcv_available <- function() {
  requireNamespace("mgcv", quietly = TRUE)
}

# qAOPweaver/qAOPweaverR/R/qker_user_data_workflow.R lines 58-192
qker_fit_user_data_candidates <- function(ke_response_data, path_structure_data) {
  validation <- qker_validate_data_intake_inputs(ke_response_data, path_structure_data)
  if (!isTRUE(validation$validation_pass)) stop("Data Intake-compatible qKER input validation failed.")
  if (!"time_point" %in% names(ke_response_data)) ke_response_data$time_point <- 0
  if (!"replicate_id" %in% names(ke_response_data)) ke_response_data$replicate_id <- ke_response_data$observation_id
  if (!"ker_label" %in% names(path_structure_data)) path_structure_data$ker_label <- path_structure_data$edge_id
  model_types <- c("linear_lm", "log_linear_lm", "quadratic_lm", "gam_candidate")
  rows <- list()
  for (edge_i in seq_len(nrow(path_structure_data))) {
    edge <- path_structure_data[edge_i, , drop = FALSE]
    upstream <- ke_response_data[ke_response_data$ke_id == edge$from_ke_id, , drop = FALSE]
    downstream <- ke_response_data[ke_response_data$ke_id == edge$to_ke_id, , drop = FALSE]
    time_pairs <- expand.grid(from_time = sort(unique(upstream$time_point)), to_time = sort(unique(downstream$time_point)))
    time_pairs <- time_pairs[time_pairs$from_time <= time_pairs$to_time, , drop = FALSE]
    for (tp_i in seq_len(nrow(time_pairs))) {
      up <- upstream[upstream$time_point == time_pairs$from_time[tp_i], , drop = FALSE]
      down <- downstream[downstream$time_point == time_pairs$to_time[tp_i], , drop = FALSE]
      keys <- intersect(c("dataset_id", "observation_id", "replicate_id", "dose"), intersect(names(up), names(down)))
      pair <- merge(up, down, by = keys, suffixes = c("_up", "_down"))
      if (!nrow(pair)) next
      x <- pair$response_value_up
      y <- pair$response_value_down
      for (mt in model_types) {
        n <- length(y)
        dat <- data.frame(x = x, y = y)
        fit <- NULL
        pred <- rep(NA_real_, n)
        eq <- NA_character_
        r2 <- NA_real_
        adj_r2 <- NA_real_
        aic <- NA_real_
        bic <- NA_real_
        deviance_explained <- NA_real_
        edf <- NA_real_
        fit_engine <- NA_character_
        model_status <- "fitted"
        notes <- "Generated by user-data-ready qKER workflow."
        if (mt == "linear_lm") {
          fit <- stats::lm(y ~ x)
          pred <- stats::fitted(fit)
          eq <- sprintf("y = %.4f + %.4f x", stats::coef(fit)[1], stats::coef(fit)[2])
          fit_engine <- "stats::lm"
        } else if (mt == "log_linear_lm") {
          fit <- stats::lm(y ~ log1p(x))
          pred <- stats::fitted(fit)
          eq <- sprintf("y = %.4f + %.4f log1p(x)", stats::coef(fit)[1], stats::coef(fit)[2])
          fit_engine <- "stats::lm"
        } else if (mt == "quadratic_lm") {
          fit <- stats::lm(y ~ x + I(x^2))
          pred <- stats::fitted(fit)
          cf <- stats::coef(fit)
          eq <- sprintf("y = %.4f + %.4f x + %.4f x^2", cf[1], cf[2], cf[3])
          fit_engine <- "stats::lm"
        } else if (!qker_v0320_mgcv_available()) {
          model_status <- "optional_candidate_unavailable"
          eq <- "GAM candidate unavailable because mgcv is not installed; no linear fallback was fitted."
          notes <- "mgcv is required for gam_candidate; qKER v0320 fails closed instead of returning a linear placeholder."
          warning("Skipping qKER gam_candidate because mgcv is not available; no linear fallback was fitted.", call. = FALSE)
        } else {
          fit <- tryCatch(mgcv::gam(y ~ s(x), data = dat), error = function(e) e)
          if (inherits(fit, "error")) {
            model_status <- "optional_candidate_unavailable_or_failed"
            eq <- "GAM candidate unavailable because mgcv::gam(y ~ s(x)) failed; no linear fallback was fitted."
            notes <- paste("mgcv::gam failed for gam_candidate:", conditionMessage(fit))
            warning(notes, call. = FALSE)
          } else {
            pred <- as.numeric(stats::predict(fit, newdata = dat))
            eq <- "GAM smooth fitted with mgcv::gam(y ~ s(x)); equation represented by the fitted smooth curve."
            fit_engine <- "mgcv::gam"
            s <- summary(fit)
            deviance_explained <- if (!is.null(s$dev.expl)) unname(s$dev.expl) else NA_real_
            edf <- if (!is.null(s$edf)) sum(s$edf) else NA_real_
          }
        }
        if (!is.null(fit) && !inherits(fit, "error")) {
          residuals <- y - pred
          sst <- sum((y - mean(y))^2)
          r2 <- if (is.finite(sst) && sst > 0) 1 - sum(residuals^2) / sst else NA_real_
          if (mt != "gam_candidate") {
            fit_summary <- suppressWarnings(summary(fit))
            adj_r2 <- suppressWarnings(fit_summary$adj.r.squared)
          }
          aic <- suppressWarnings(stats::AIC(fit))
          bic <- suppressWarnings(stats::BIC(fit))
        }
        rmse <- if (all(is.na(pred))) NA_real_ else sqrt(mean((y - pred)^2))
        mae <- if (all(is.na(pred))) NA_real_ else mean(abs(y - pred))
        rows[[length(rows) + 1]] <- data.frame(
          run_id = paste0(unique(ke_response_data$dataset_id), "_user_data_ready_qker_run"),
          dataset_id = unique(ke_response_data$dataset_id),
          edge_id = edge$edge_id,
          ker_label = edge$ker_label,
          from_ke_id = edge$from_ke_id,
          to_ke_id = edge$to_ke_id,
          from_ke_label = if ("from_ke_label" %in% names(edge)) edge$from_ke_label else edge$from_ke_id,
          to_ke_label = if ("to_ke_label" %in% names(edge)) edge$to_ke_label else edge$to_ke_id,
          from_time = time_pairs$from_time[tp_i],
          to_time = time_pairs$to_time[tp_i],
          pair_label = paste(time_pairs$from_time[tp_i], "->", time_pairs$to_time[tp_i]),
          model_id = paste0(edge$edge_id, "_", mt, "_t", gsub("[.]", "p", time_pairs$from_time[tp_i]), "_t", gsub("[.]", "p", time_pairs$to_time[tp_i])),
          model_type = mt,
          model_family = sub("_lm$", "", mt),
          n = n,
          r_squared = round(r2, 6),
          adjusted_r_squared = round(adj_r2, 6),
          rmse = rmse,
          mae = mae,
          aic = aic,
          bic = bic,
          deviance_explained = deviance_explained,
          edf = edf,
          fit_engine = fit_engine,
          p_value = NA,
          model_equation = eq,
          equation_display = eq,
          candidate_rank = NA_integer_,
          candidate_scope = "all_valid_time_pair_model_candidates",
          candidate_source = "Data Intake-compatible user data",
          model_status = model_status,
          notes = notes,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  candidates <- do.call(rbind, rows)
  candidates <- candidates[order(candidates$edge_id, candidates$from_time, candidates$to_time, candidates$rmse, na.last = TRUE), , drop = FALSE]
  candidates$candidate_rank <- ave(candidates$rmse, candidates$edge_id, candidates$from_time, candidates$to_time, FUN = function(x) {
    out <- rep(NA_integer_, length(x))
    ok <- is.finite(x)
    out[ok] <- rank(x[ok], ties.method = "first")
    out
  })
  candidates
}

# qAOPweaver/qAOPweaverR/R/ode_contract.R lines 60-67
ode_empty_diag <- function() {
  data.frame(
    code = character(),
    field = character(),
    message = character(),
    stringsAsFactors = FALSE
  )
}

# qAOPweaver/qAOPweaverR/R/ode_contract.R lines 69-76
ode_diag <- function(code, field, message) {
  data.frame(
    code = as.character(code),
    field = as.character(field),
    message = as.character(message),
    stringsAsFactors = FALSE
  )
}

# qAOPweaver/qAOPweaverR/R/ode_contract.R lines 78-82
ode_bind_diag <- function(items) {
  items <- Filter(function(x) is.data.frame(x) && nrow(x), items)
  if (!length(items)) return(ode_empty_diag())
  do.call(rbind, items)
}

# qAOPweaver/qAOPweaverR/R/ode_model_spec.R lines 18-20
ode_model_allowed_math_functions <- function() {
  c("exp", "log", "log10", "sqrt", "sin", "cos", "tan", "tanh", "pmin", "pmax", "abs", "min", "max")
}

# qAOPweaver/qAOPweaverR/R/ode_simulation.R lines 10-12
ode_sim_empty_diag <- function() {
  ode_empty_diag()
}

# qAOPweaver/qAOPweaverR/R/ode_simulation.R lines 14-16
ode_sim_diag <- function(code, field, message) {
  ode_diag(code, field, message)
}

# qAOPweaver/qAOPweaverR/R/ode_simulation.R lines 18-20
ode_sim_bind_diag <- function(items) {
  ode_bind_diag(items)
}

# qAOPweaver/qAOPweaverR/R/ode_simulation.R lines 262-283
ode_sim_eval_derivatives <- function(state_values, time, equations, parameters, forcing_schedule = data.frame()) {
  allowed_functions <- stats::setNames(lapply(ode_model_allowed_math_functions(), get, envir = baseenv()), ode_model_allowed_math_functions())
  allowed_operators <- stats::setNames(lapply(c("+", "-", "*", "/", "^", "("), get, envir = baseenv()), c("+", "-", "*", "/", "^", "("))
  out <- stats::setNames(rep(NA_real_, length(state_values)), names(state_values))
  forcing <- ode_sim_forcing_values_at(forcing_schedule, time)
  if (nrow(forcing$errors)) return(list(values = out, errors = forcing$errors))
  for (state_id in names(state_values)) {
    expr <- equations[[state_id]]
    env <- list2env(as.list(c(state_values, parameters, forcing$values, time = time)), parent = emptyenv())
    for (fn in names(allowed_functions)) assign(fn, allowed_functions[[fn]], envir = env)
    for (op in names(allowed_operators)) assign(op, allowed_operators[[op]], envir = env)
    value <- tryCatch(eval(parse(text = expr), envir = env), error = function(e) e)
    if (inherits(value, "error") || length(value) != 1L) {
      return(list(values = out, errors = ode_sim_diag("solver_failure", paste0("equations.", state_id), if (inherits(value, "error")) conditionMessage(value) else "Derivative expression did not return a scalar.")))
    }
    out[[state_id]] <- suppressWarnings(as.numeric(value))
  }
  if (any(!is.finite(out))) {
    return(list(values = out, errors = ode_sim_diag("non_finite_derivative", "equations", "One or more derivative values were non-finite.")))
  }
  list(values = out, errors = ode_sim_empty_diag())
}

# Adapted from qAOPweaver/qAOPweaverR/R/ode_simulation.R lines 285-320 (RK4 only)
ode_sim_integrate <- function(model_object, parameters, initial_values, times, method, forcing_schedule = data.frame()) {
  stopifnot(identical(method, "rk4"))
  eq <- model_object$equations
  active <- eq[is.na(eq$is_active) | eq$is_active, , drop = FALSE]
  equations <- stats::setNames(as.character(active$expression), as.character(active$target_state))
  states <- names(initial_values)
  y <- initial_values
  rows <- list(data.frame(time = times[[1]], state_id = states, simulated_value = as.numeric(y), stringsAsFactors = FALSE))
  errors <- list()
  deriv <- function(vals, t) ode_sim_eval_derivatives(vals, t, equations, parameters, forcing_schedule = forcing_schedule)
  for (i in seq_len(length(times) - 1L)) {
    t0 <- times[[i]]
    t1 <- times[[i + 1L]]
    h <- t1 - t0
    k1 <- deriv(y, t0)
    if (nrow(k1$errors)) { errors <- c(errors, list(k1$errors)); break }
    k2 <- deriv(y + h * k1$values / 2, t0 + h / 2)
    if (nrow(k2$errors)) { errors <- c(errors, list(k2$errors)); break }
    k3 <- deriv(y + h * k2$values / 2, t0 + h / 2)
    if (nrow(k3$errors)) { errors <- c(errors, list(k3$errors)); break }
    k4 <- deriv(y + h * k3$values, t1)
    if (nrow(k4$errors)) { errors <- c(errors, list(k4$errors)); break }
    y <- y + h * (k1$values + 2 * k2$values + 2 * k3$values + k4$values) / 6
    if (any(!is.finite(y))) {
      errors <- c(errors, list(ode_sim_diag("non_finite_simulated_state", "simulation_results", "One or more simulated state values were non-finite.")))
      break
    }
    rows[[length(rows) + 1L]] <- data.frame(time = t1, state_id = states, simulated_value = as.numeric(y), stringsAsFactors = FALSE)
  }
  list(results = do.call(rbind, rows), errors = ode_sim_bind_diag(errors))
}

# qAOPweaver/qAOPweaverR/R/ode_fitting.R lines 231-265
ode_fit_align_predictions <- function(simulation_results, observed_data) {
  sim <- simulation_results
  if (!nrow(sim)) {
    return(list(alignment = data.frame(n_observations = nrow(observed_data), n_matched = 0L, n_unmatched = nrow(observed_data), stringsAsFactors = FALSE), predictions = data.frame(), residuals = data.frame(), errors = ode_diag("prediction_alignment_failed", "fitted_predictions", "No simulation results were available for prediction alignment.")))
  }
  sim$key <- paste(sim$time, sim$state_id, sim$scenario_id, sep = "::")
  obs <- observed_data
  obs$key <- paste(obs$time, obs$state_id, obs$scenario_id, sep = "::")
  pred <- sim[match(obs$key, sim$key), , drop = FALSE]
  missing <- is.na(pred$simulated_value)
  if (any(missing)) {
    alignment <- data.frame(n_observations = nrow(obs), n_matched = sum(!missing), n_unmatched = sum(missing), stringsAsFactors = FALSE)
    return(list(alignment = alignment, predictions = data.frame(), residuals = data.frame(), errors = ode_diag("prediction_alignment_failed", "observed_data", "One or more observed time/state/scenario keys were not matched to predictions.")))
  }
  predictions <- data.frame(
    observation_id = obs$observation_id,
    time = obs$time,
    state_id = obs$state_id,
    scenario_id = obs$scenario_id,
    replicate_id = obs$replicate_id,
    observed_value = obs$observed_value,
    predicted_value = pred$simulated_value,
    stringsAsFactors = FALSE
  )
  residuals <- predictions
  residuals$residual <- residuals$observed_value - residuals$predicted_value
  residuals$absolute_residual <- abs(residuals$residual)
  residuals$squared_residual <- residuals$residual^2
  list(
    alignment = data.frame(n_observations = nrow(obs), n_matched = nrow(obs), n_unmatched = 0L, alignment_key = "time;state_id;scenario_id", stringsAsFactors = FALSE),
    predictions = predictions,
    residuals = residuals,
    errors = ode_empty_diag()
  )
}

# qAOPweaver/qAOPweaverR/R/ode_fitting.R lines 267-285
ode_fit_metrics <- function(residuals, n_fitted_parameters, objective_value) {
  if (!nrow(residuals)) {
    return(data.frame(RMSE = NA_real_, MAE = NA_real_, bias = NA_real_, objective_value = objective_value, n_observations = 0L, n_fitted_parameters = n_fitted_parameters, AIC = NA_real_, BIC = NA_real_, metric_scope = "training_fit_only", stringsAsFactors = FALSE))
  }
  rss <- sum(residuals$squared_residual)
  n <- nrow(residuals)
  data.frame(
    RMSE = sqrt(mean(residuals$squared_residual)),
    MAE = mean(residuals$absolute_residual),
    bias = mean(residuals$residual),
    objective_value = objective_value %||% rss,
    n_observations = n,
    n_fitted_parameters = n_fitted_parameters,
    AIC = NA_real_,
    BIC = NA_real_,
    metric_scope = "training_fit_only",
    stringsAsFactors = FALSE
  )
}
