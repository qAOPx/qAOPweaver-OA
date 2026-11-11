# Functions copied from qAOPweaver/qAOPweaverR/R/prepared_analysis_v1.R
# qAOPweaver v0.29.0 revision 0143dd56dd21a5268b17ecbe8383647b1c6323c9.
# See the parent repository's provenance and licensing documentation.

pa_fail <- function(...) stop(paste0(...), call. = FALSE)

pa_value <- function(x, default) if (is.null(x)) default else x

pa_condition_cluster_inference <- function(x, y, weights, coefficients, condition_ids) {
  groups <- split(seq_along(y), as.character(condition_ids))
  n_groups <- length(groups)
  n_parameters <- ncol(x)
  unavailable <- function(estimate, reason) list(
    status = "not_available", method = "condition_cluster_sandwich_t",
    estimate = as.numeric(estimate), independent_conditions = n_groups,
    standard_error = NULL, p_value = NULL, confidence_interval_95 = NULL,
    reason = reason)
  result <- setNames(vector("list", n_parameters), colnames(x))
  if (n_groups <= n_parameters) {
    for (term in colnames(x)) result[[term]] <- unavailable(
      coefficients[[term]], "Too few independent conditions for this component regression")
    return(result)
  }
  residual <- as.numeric(y - x %*% coefficients)
  bread <- tryCatch(solve(crossprod(x, weights * x)), error = function(e) NULL)
  if (is.null(bread) || any(!is.finite(bread))) {
    for (term in colnames(x)) result[[term]] <- unavailable(
      coefficients[[term]], "Condition-cluster covariance is singular")
    return(result)
  }
  meat <- matrix(0, n_parameters, n_parameters)
  for (indices in groups) {
    score <- as.numeric(crossprod(x[indices, , drop = FALSE],
      weights[indices] * residual[indices]))
    meat <- meat + tcrossprod(score)
  }
  covariance <- (n_groups / (n_groups - 1)) * bread %*% meat %*% bread
  degrees_freedom <- n_groups - 1L
  for (term in colnames(x)) {
    index <- match(term, colnames(x))
    estimate <- as.numeric(coefficients[[term]])
    variance <- covariance[index, index]
    if (!is.finite(variance) || variance <= 0) {
      result[[term]] <- unavailable(estimate,
        "Condition-cluster variance is zero or unavailable")
      next
    }
    se <- sqrt(variance)
    critical <- stats::qt(.975, df = degrees_freedom)
    result[[term]] <- list(
      status = "exploratory", method = "condition_cluster_sandwich_t",
      estimate = estimate, independent_conditions = n_groups,
      standard_error = se,
      p_value = 2 * stats::pt(-abs(estimate / se), df = degrees_freedom),
      confidence_interval_95 = as.list(c(estimate - critical * se,
        estimate + critical * se)),
      degrees_freedom = degrees_freedom,
      note = paste("Unadjusted exploratory condition-cluster inference under",
        "the selected Step 2 pairing; constructed rows are not independent",
        "biological replicates"))
  }
  result
}

pa_request <- function(request, model_type) {
  if (!identical(request$schema_version, "prepared_analysis_request_v1") ||
      !identical(request$model_type, model_type))
    pa_fail("Expected prepared_analysis_request_v1 for ", model_type)
  view <- request$prepared_view
  if (!identical(view$manifest$schema_version, "step2_prepared_view_v1") ||
      length(view$manifest$generation) != 1L ||
      length(view$manifest$view_digest) != 1L ||
      !nzchar(as.character(view$manifest$view_digest)))
    pa_fail("A versioned Step 2 generation and digest are required")
  graph <- bs_graph(request$network)
  nodes <- setNames(graph$nodes, graph$ids)
  for (edge in graph$edges) {
    from <- nodes[[edge$from]]; to <- nodes[[edge$to]]
    ft <- pa_value(from$fit_time, NULL); tt <- pa_value(to$fit_time, NULL)
    if (!is.null(ft) || !is.null(tt)) {
      if (is.null(ft) || is.null(tt) || length(ft) != 1L || length(tt) != 1L ||
          !is.numeric(ft) || !is.numeric(tt) || !is.finite(ft) || !is.finite(tt))
        pa_fail("Every temporal edge requires finite numeric times at both nodes")
      from_id <- as.character(pa_value(from$node_id, from$id))
      to_id <- as.character(pa_value(to$node_id, to$id))
      if (identical(from_id, to_id) && ft >= tt)
        pa_fail("A self-memory edge requires strictly increasing time")
      if (ft > tt) pa_fail("An AOP upstream node cannot have a later time than its downstream node")
    }
  }
  cases <- bs_frame(view$cases)
  needed <- c("case_id", "condition_id", "case_weight", graph$ids)
  if (!nrow(cases) || !all(needed %in% names(cases)))
    pa_fail("Prepared cases require ", paste(needed, collapse = ", "))
  if (nrow(cases) > 100000L || anyNA(cases$case_id) || anyDuplicated(cases$case_id) ||
      anyNA(cases$condition_id) || any(!nzchar(as.character(cases$condition_id))))
    pa_fail("Prepared case/condition identities must be unique and complete")
  for (node in graph$ids) {
    if (!is.numeric(cases[[node]]) || any(!is.finite(cases[[node]])))
      pa_fail("Prepared node ", node, " needs finite numeric values in every case")
  }
  weights <- cases$case_weight
  if (!is.numeric(weights) || any(!is.finite(weights)) || any(weights <= 0))
    pa_fail("Prepared case weights must be finite and positive")
  for (key in c("run_id", "realization_id")) if (key %in% names(cases) &&
      length(unique(cases[[key]])) > 1L)
    pa_fail("Fit one saved ", key, " at a time; never stack pairing realizations")
  kind <- as.character(pa_value(view$manifest$view_kind, ""))
  if (kind %in% c("constructed_cases", "hypothetical_cases", "cartesian_cases")) {
    mass <- tapply(weights, as.character(cases$condition_id), sum)
    if (any(abs(mass - 1) > 1e-8))
      pa_fail("Each constructed condition must have total case weight one")
  }
  list(view = view, graph = graph, cases = cases, weights = weights,
    view_kind = kind)
}

run_prepared_psem_v1 <- function(request) {
  input <- pa_request(request, "psem")
  config <- pa_value(request$config, list())
  family <- pa_value(config$component_family, "weighted_linear")
  if (!family %in% c("weighted_linear", "linear"))
    pa_fail("Prepared PSEM currently supports weighted linear components only; ",
      "mixed and GLM families are unavailable for this view")
  graph <- input$graph; cases <- input$cases; weights <- input$weights
  equations <- list()
  for (node in graph$order) {
    parents <- graph$parents[[node]]
    if (!length(parents)) next
    x <- cbind(`(Intercept)` = 1, as.matrix(cases[, parents, drop = FALSE]))
    colnames(x) <- c("(Intercept)", parents)
    y <- cases[[node]]
    fit <- stats::lm.wfit(x = x, y = y, w = weights)
    if (fit$rank < ncol(x) || any(!is.finite(fit$coefficients)) ||
        nrow(cases) <= ncol(x))
      pa_fail("Component regression is singular or lacks residual degrees of freedom for ", node)
    predicted <- as.numeric(x %*% fit$coefficients)
    rss <- sum(weights * (y - predicted)^2)
    center <- sum(weights * y) / sum(weights)
    tss <- sum(weights * (y - center)^2)
    equations[[node]] <- list(
      node_id = node, parent_ids = as.list(parents),
      coefficients = as.list(stats::setNames(as.numeric(fit$coefficients), colnames(x))),
      path_inference = pa_condition_cluster_inference(x, y, weights,
        fit$coefficients, cases$condition_id),
      n_cases = nrow(cases), n_conditions = length(unique(cases$condition_id)),
      weight_sum = sum(weights), weighted_rmse = sqrt(rss / sum(weights)),
      weighted_r_squared = if (tss > 0) 1 - rss / tss else NULL,
      fitted = as.list(predicted), residuals = as.list(y - predicted)
    )
  }
  if (!length(equations)) pa_fail("PSEM requires at least one endogenous node")
  constructed <- input$view_kind %in% c("constructed_cases", "hypothetical_cases", "cartesian_cases")
  list(
    schema_version = "prepared_psem_model_v1", status = "success",
    view_digest = as.character(input$view$manifest$view_digest),
    generation = input$view$manifest$generation,
    view_kind = input$view_kind, network = request$network,
    equations = equations,
    inference_scope = if (constructed)
      paste("Exploratory constructed cases; condition-cluster uncertainty and",
        "p-values are conditional on the selected pairing, not independent-row",
        "or biological-replicate significance")
    else paste("Exploratory condition-cluster uncertainty on selected observed or",
      "condition-level cases; small numbers of independent conditions limit inference"),
    provenance = list(case_ids = as.list(as.character(cases$case_id)),
      condition_ids = as.list(as.character(cases$condition_id)),
      case_weights = as.list(weights))
  )
}

predict_prepared_psem_v1 <- function(request) {
  model <- request$model
  if (!identical(model$schema_version, "prepared_psem_model_v1"))
    pa_fail("Expected a frozen prepared PSEM model")
  graph <- bs_graph(model$network)
  if (!setequal(names(model$equations), names(graph$parents)[lengths(graph$parents) > 0L]))
    pa_fail("Frozen PSEM equations do not match its graph")
  mode <- pa_value(request$mode, "one_step")
  if (!identical(mode, "one_step")) pa_fail("The case uses one-step prediction")
  cases <- bs_frame(request$cases)
  if (!nrow(cases)) pa_fail("Prediction requires cases")
  if (!"case_id" %in% names(cases) || anyNA(cases$case_id) || anyDuplicated(cases$case_id))
    pa_fail("Prediction case_id must be unique and complete")
  targets <- bs_vector(pa_value(request$targets, as.list(names(model$equations))))
  if (!length(targets) || any(!targets %in% names(model$equations)) || anyDuplicated(targets))
    pa_fail("Prediction targets must be unique modeled endogenous nodes")
  required <- unique(unlist(graph$parents[targets], use.names = FALSE))
  if (!all(required %in% names(cases)))
    pa_fail("Prediction is missing required observed input columns: ",
      paste(setdiff(required, names(cases)), collapse = ", "))
  for (node in required) if (!is.numeric(cases[[node]]) || any(!is.finite(cases[[node]])))
    pa_fail("Prediction input ", node, " requires finite numeric observations")
  predictions <- list()
  for (node in graph$order[graph$order %in% targets]) {
    equation <- model$equations[[node]]
    if (is.null(equation)) next
    parents <- bs_vector(equation$parent_ids)
    coefficients <- unlist(equation$coefficients, use.names = TRUE)
    if (!setequal(names(coefficients), c("(Intercept)", parents)) ||
        any(!is.finite(coefficients)))
      pa_fail("Frozen equation coefficients are invalid for ", node)
    if (!all(parents %in% names(cases)))
      pa_fail("Prediction lacks a required parent value for ", node)
    estimate <- rep(as.numeric(coefficients[["(Intercept)"]]), nrow(cases))
    for (parent in parents) estimate <- estimate +
      as.numeric(coefficients[[parent]]) * cases[[parent]]
    predictions[[node]] <- estimate
  }
  list(schema_version = "prepared_psem_prediction_v1", status = "success",
    model_view_digest = model$view_digest, mode = mode,
    case_ids = as.list(as.character(cases$case_id)),
    predictions = setNames(lapply(targets, function(node)
      as.list(as.numeric(predictions[[node]]))), targets),
    interpretation = "Uses measured direct-parent inputs; target outcomes are never inputs")
}

score_prepared_psem_v1 <- function(request) {
  model <- request$model
  target <- as.character(request$target)
  if (length(target) != 1L || !target %in% names(model$equations))
    pa_fail("Scoring requires one modeled PSEM target")
  cases <- bs_frame(request$cases)
  if (!nrow(cases) || !target %in% names(cases) ||
      !is.numeric(cases[[target]]) || any(!is.finite(cases[[target]])))
    pa_fail("External scoring requires finite observed target values")
  weights <- if ("case_weight" %in% names(cases)) cases$case_weight else rep(1, nrow(cases))
  if (!is.numeric(weights) || any(!is.finite(weights)) || any(weights <= 0))
    pa_fail("Scoring weights must be finite and positive")
  predicted <- predict_prepared_psem_v1(list(model = model, cases = cases,
    mode = pa_value(request$mode, "one_step"), targets = target))
  values <- as.numeric(unlist(predicted$predictions[[target]], use.names = FALSE))
  residuals <- cases[[target]] - values
  list(schema_version = "prepared_psem_external_score_v1", status = "success",
    model_view_digest = model$view_digest, target = target, mode = predicted$mode,
    scored_cases = nrow(cases), weight_sum = sum(weights),
    rmse = sqrt(sum(weights * residuals^2) / sum(weights)),
    mae = sum(weights * abs(residuals)) / sum(weights),
    bias = sum(weights * residuals) / sum(weights),
    predictions = lapply(seq_len(nrow(cases)), function(i) list(
      case_id = as.character(cases$case_id[i]), target = target,
      observed = cases[[target]][i], predicted = values[i],
      residual = residuals[i], case_weight = weights[i])),
    interpretation = "Frozen component equations; no external refit, retuning or target-outcome input")
}

pa_weighted_quantile <- function(x, w, p) {
  values <- sort(unique(x))
  if (length(values) < 2L) pa_fail("Quantile states need at least two distinct response values")
  mass <- vapply(values, function(value) sum(w[x == value]), numeric(1))
  index <- which(cumsum(mass) / sum(mass) >= p)[1]
  index <- min(index, length(values) - 1L)
  values[index] + (values[index + 1L] - values[index]) / 2
}

pa_state_design <- function(input, config) {
  graph <- input$graph; cases <- input$cases; weights <- input$weights
  nodes <- setNames(graph$nodes, graph$ids)
  base <- vapply(graph$ids, function(n) as.character(pa_value(nodes[[n]]$node_id, n)), character(1))
  design <- pa_value(config$state_design, list())
  states <- list(); cuts <- list(); values <- list(); intervals <- list(); info <- list()
  for (base_id in unique(base)) {
    members <- graph$ids[base == base_id]
    cfg <- pa_value(design[[base_id]], pa_value(design[[members[1]]], list()))
    method <- pa_value(cfg$method, "quantile")
    n_states <- bs_number(cfg$n_states, 3, "n_states", 2, 10, TRUE)
    x <- unlist(lapply(members, function(n) cases[[n]]), use.names = FALSE)
    w <- rep(weights, times = length(members))
    fixed <- pa_value(cfg$fixed_cutpoints, NULL)
    if (!is.null(fixed)) {
      if (!method %in% c("quantile", "equal_width"))
        pa_fail("Fixed state cutpoints require quantile or equal_width origin")
      cut <- unlist(fixed, use.names = FALSE)
      if (!is.numeric(cut) || length(cut) != n_states - 1L ||
          any(!is.finite(cut)) || any(diff(cut) <= 0))
        pa_fail("Frozen state cutpoints must be finite, strictly increasing and match n_states")
    } else if (method == "equal_width") {
      cut <- seq(min(x), max(x), length.out = n_states + 1L)[-c(1L, n_states + 1L)]
    } else if (method == "quantile") {
      cut <- vapply(seq_len(n_states - 1L) / n_states,
        function(p) pa_weighted_quantile(x, w, p), numeric(1))
    } else pa_fail("Prepared-case state method is unsupported: ", method)
    if (anyDuplicated(cut)) pa_fail("Selected data cannot realize distinct state cutpoints for ", base_id)
    for (node in members) {
      labels <- paste0("state_", seq_len(n_states))
      states[[node]] <- as.list(labels); cuts[[node]] <- cut
      values[[node]] <- bs_assign(cases[[node]], cut)
      bounds <- c(-Inf, cut, Inf)
      for (i in seq_len(n_states)) intervals[[length(intervals) + 1L]] <- list(
        node_id = node, state_id = labels[i],
        lower = if (is.finite(bounds[i])) bounds[i] else NULL,
        upper = if (is.finite(bounds[i + 1L])) bounds[i + 1L] else NULL,
        lower_inclusive = is.finite(bounds[i]), upper_inclusive = FALSE,
        weighted_support = sum(weights[values[[node]] == labels[i]]),
        case_count = sum(values[[node]] == labels[i]))
      info[[length(info) + 1L]] <- list(node_id = node, base_node_id = base_id,
        method = method, n_states = n_states, cutpoints = as.list(cut),
        state_basis = if (!is.null(fixed))
          "Frozen reference cutpoints shared across pairing realizations" else
          "Shared across all selected times for this KE")
    }
  }
  list(states = states, cuts = cuts, values = values,
    discretization = list(intervals = intervals, nodes = info))
}

pa_cpts <- function(input, design, config) {
  graph <- input$graph; values <- design$values; weights <- input$weights
  alpha <- bs_number(config$smoothing_alpha, 1, "smoothing_alpha", 0)
  if (length(config$root_priors)) pa_fail("The case estimates empirical root probabilities")
  cpts <- list(); prior_records <- list(); warnings <- character()
  for (node in graph$ids) {
    parents <- graph$parents[[node]]; labels <- bs_vector(design$states[[node]])
    rows <- lapply(bs_assignments(parents, design$states), function(assignment) {
      take <- rep(TRUE, length(weights))
      for (parent in parents) take <- take & values[[parent]] == assignment[[parent]]
      counts <- vapply(labels, function(label) sum(weights[take & values[[node]] == label]), numeric(1))
      names(counts) <- labels
      n_cases <- sum(take)
      if (!length(parents)) {
        method <- "empirical"
        if (sum(counts) + length(labels) * alpha <= 0)
          pa_fail("Empty root support with zero smoothing for ", node)
        probability <- (counts + alpha) / (sum(counts) + length(labels) * alpha)
        prior_records[[node]] <<- list(method = method,
          probabilities = setNames(as.list(probability), labels),
          weighted_counts = as.list(counts))
      } else {
        if (sum(counts) + length(labels) * alpha <= 0)
          pa_fail("Empty parent configuration with zero smoothing for ", node)
        probability <- (counts + alpha) / (sum(counts) + length(labels) * alpha)
      }
      if (n_cases == 0) warnings <<- c(warnings,
        paste0("No selected case for one parent configuration of ", node))
      list(parent_assignment = assignment,
        probabilities = setNames(as.list(as.numeric(probability)), labels),
        counts = as.list(counts), weighted_support = sum(counts),
        case_count = n_cases)
    })
    cpts[[node]] <- list(node_id = node, parent_ids = as.list(parents), rows = rows,
      source = if (length(parents)) "prepared_weighted_cases" else "root_prior",
      counting_unit = "Selected Step 2 case weight, not independent biological replicates",
      smoothing_alpha = alpha,
      smoothing_rule = "(weighted_state_support + alpha) / (weighted_family_support + states * alpha)")
  }
  list(cpts = cpts, priors = prior_records, warnings = as.list(unique(warnings)), alpha = alpha)
}

pa_regression_cpts <- function(input, design, counted, config, mode = "regression") {
  graph <- input$graph; cases <- input$cases; weights <- input$weights
  bayesian <- identical(mode, "bayesian_regression")
  simulation <- pa_value(pa_value(config$regression, list())$simulation, list())
  count <- bs_number(simulation$n, 5000, "simulation n", 100, 100000, TRUE)
  seed <- bs_number(simulation$seed, 42, "simulation seed", 0, .Machine$integer.max, TRUE)
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit(if (is.null(old_seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  } else assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
  set.seed(seed)

  cpts <- counted$cpts; details <- list()
  for (node in graph$ids) {
    parents <- graph$parents[[node]]
    if (!length(parents)) next
    parent <- parents[[1]]
    x <- cases[[parent]]; y <- cases[[node]]
    matrix <- cbind(`(Intercept)` = 1, parent = x)
    fit <- stats::lm.wfit(matrix, y, weights)
    if (nrow(cases) < 3L || fit$rank != 2L || any(!is.finite(fit$coefficients)))
      pa_fail("Prepared regression CPT requires predictor variation and a full-rank weighted fit for ", node)
    residual_sd <- sqrt(sum(weights * fit$residuals^2) / sum(weights))
    if (!is.finite(residual_sd)) pa_fail("Prepared regression CPT has nonfinite residual scatter for ", node)
    if (bayesian) {
      effective_weight <- sum(weights)
      degrees_freedom <- effective_weight - 2
      if (length(unique(cases$condition_id)) < 3L || degrees_freedom <= 0)
        pa_fail("Posterior CPT requires more than two effective weight units and at least three independent conditions for ", node)
      residual_sum_squares <- sum(weights * fit$residuals^2)
      if (!is.finite(residual_sum_squares) ||
          residual_sum_squares <= .Machine$double.eps * max(1, sum(weights * y^2)))
        pa_fail("Posterior CPT requires positive residual variation for ", node)
      coefficient_covariance <- tryCatch(solve(crossprod(matrix, weights * matrix)),
        error = function(error) NULL)
      if (is.null(coefficient_covariance) || any(!is.finite(coefficient_covariance)))
        pa_fail("Posterior CPT weighted coefficient covariance is singular for ", node)
      coefficient_cholesky <- tryCatch(chol(coefficient_covariance),
        error = function(error) NULL)
      if (is.null(coefficient_cholesky))
        pa_fail("Posterior CPT weighted coefficient covariance is not positive definite for ", node)
      parameter_draws <- vector("list", length(design$states[[parent]]))
      variance_draws <- vector("list", length(design$states[[parent]]))
    }

    parent_states <- bs_vector(design$states[[parent]])
    child_states <- bs_vector(design$states[[node]])
    allocation <- rep(count %/% length(parent_states), length(parent_states))
    if (count %% length(parent_states))
      allocation[seq_len(count %% length(parent_states))] <-
        allocation[seq_len(count %% length(parent_states))] + 1L
    rows <- vector("list", length(parent_states))
    for (i in seq_along(parent_states)) {
      selected <- design$values[[parent]] == parent_states[i]
      if (!any(selected))
        pa_fail("Regression CPT parent state has no selected Step 2 case for ", parent,
          " state ", parent_states[i])
      selected_x <- x[selected]
      draw_x <- selected_x[sample.int(length(selected_x), size = allocation[i],
        replace = TRUE, prob = weights[selected])]
      if (bayesian) {
        variance <- residual_sum_squares / stats::rchisq(allocation[i], df = degrees_freedom)
        if (any(!is.finite(variance)) || any(variance <= 0))
          pa_fail("Posterior CPT variance draw is nonfinite for ", node)
        normals <- matrix(stats::rnorm(2L * allocation[i]), ncol = 2L)
        coefficients <- sweep(normals %*% coefficient_cholesky, 1L,
          sqrt(variance), `*`) +
          matrix(rep(as.numeric(fit$coefficients), each = allocation[i]), ncol = 2L)
        draw_y <- coefficients[, 1] + coefficients[, 2] * draw_x +
          stats::rnorm(allocation[i], sd = sqrt(variance))
        parameter_draws[[i]] <- coefficients
        variance_draws[[i]] <- variance
      } else {
        draw_y <- fit$coefficients[1] + fit$coefficients[2] * draw_x +
          stats::rnorm(allocation[i], sd = residual_sd)
      }
      if (any(!is.finite(draw_y))) pa_fail("Nonfinite regression CPT simulation for ", node)
      simulated <- tabulate(match(bs_assign(draw_y, design$cuts[[node]]), child_states),
        nbins = length(child_states))
      names(simulated) <- child_states
      probability <- (simulated + counted$alpha) /
        (sum(simulated) + length(child_states) * counted$alpha)
      source <- counted$cpts[[node]]$rows[[i]]
      rows[[i]] <- list(parent_assignment = source$parent_assignment,
        probabilities = setNames(as.list(as.numeric(probability)), child_states),
        counts = as.list(simulated),
        simulated_counts = as.list(simulated),
        simulated_draw_count = allocation[i],
        source_weighted_counts = source$counts,
        source_weighted_support = source$weighted_support,
        source_case_count = source$case_count,
        source_condition_count = length(unique(cases$condition_id[selected])))
    }
    cpts[[node]]$rows <- rows
    cpts[[node]]$source <- if (bayesian)
      "prepared_weighted_bayesian_regression_posterior_predictive" else
      "prepared_weighted_regression_residual_scatter"
    cpts[[node]]$counting_unit <- "Seeded simulated draw; measured Step 2 support is reported separately"
    cpts[[node]]$smoothing_rule <-
      "(simulated_state_count + alpha) / (simulated_parent_state_draws + states * alpha)"
    edge <- list(from = parent, to = node,
      family = "weighted_linear", coefficients = as.list(fit$coefficients),
      source_case_count = nrow(cases), source_condition_count =
        length(unique(cases$condition_id)), source_weight_sum = sum(weights),
      simulated_draw_count = count, seed = seed, residual_sd = residual_sd,
      residual_sd_rule = "Weighted root mean square of selected-case residuals",
      parent_draw_rule = "Weighted resampling of observed selected Step 2 parent values within each parent state",
      child_draw_rule = if (bayesian)
        "Joint posterior coefficient and variance draw, then Normal predictive residual scatter" else
        "Fitted weighted linear mean plus Normal(0, residual_sd) scatter",
      uncertainty_scope = if (bayesian)
        "Weighted Gaussian pseudo-posterior parameter and predictive residual uncertainty; no pairing, model-family, or external-data uncertainty" else
        "Residual scatter only; no fitted-mean, source-resampling, or posterior uncertainty")
    if (bayesian) {
      all_coefficients <- do.call(rbind, parameter_draws)
      all_variances <- unlist(variance_draws, use.names = FALSE)
      interval <- function(x) as.list(as.numeric(stats::quantile(x, c(.025, .975), names = FALSE)))
      edge$posterior <- list(
        prior = "reference_1_over_sigma_squared",
        likelihood = "weighted_gaussian_power_likelihood",
        effective_weight_sum = effective_weight,
        independent_conditions = length(unique(cases$condition_id)),
        degrees_freedom = degrees_freedom,
        weighted_residual_sum_squares = residual_sum_squares,
        coefficient_intervals = list(intercept = interval(all_coefficients[, 1]),
          parent = interval(all_coefficients[, 2])),
        variance_interval = interval(all_variances),
        draw_rule = "sigma_squared = weighted_SSE / chi_square(df); beta | sigma_squared is joint Normal with covariance sigma_squared * inverse(X'WX); child response has Normal residual scatter")
    }
    details[[length(details) + 1L]] <- edge
  }
  list(cpts = cpts, regression = list(status = "active", edges = details,
    total_draws_per_edge = count, seed = seed))
}

run_prepared_bn_v1 <- function(request) {
  input <- pa_request(request, "bn")
  config <- pa_value(request$config, list())
  mode <- pa_value(config$cpt_mode, "observed_counts")
  if (!mode %in% c("observed_counts", "regression", "bayesian_regression"))
    pa_fail("Prepared weighted BN CPT mode is unsupported: ", mode)
  if (mode %in% c("regression", "bayesian_regression") && any(lengths(input$graph$parents) > 1L))
    pa_fail("Regression CPT mode supports one parent per child; use observed_counts for convergent nodes")
  design <- pa_state_design(input, config)
  bs_limits(input$graph, design$states)
  counted <- pa_cpts(input, design, config)
  regression <- NULL
  if (mode %in% c("regression", "bayesian_regression")) {
    simulated <- pa_regression_cpts(input, design, counted, config, mode)
    counted$cpts <- simulated$cpts; regression <- simulated$regression
  }
  model <- list(schema_version = "bn_static_cases_model_v1",
    prepared_schema_version = "prepared_bn_model_v1", status = "success",
    view_digest = as.character(input$view$manifest$view_digest),
    generation = input$view$manifest$generation,
    view_kind = input$view_kind,
    bn = list(nodes = input$graph$nodes, edges = input$graph$edges,
      node_states = design$states, cpts = counted$cpts),
    discretization = design$discretization, priors = counted$priors,
    warnings = counted$warnings,
    config = list(cpt_mode = mode, smoothing_alpha = counted$alpha,
      state_design = config$state_design),
    provenance = list(case_ids = as.list(as.character(input$cases$case_id)),
      condition_ids = as.list(as.character(input$cases$condition_id)),
      case_weights = as.list(input$weights),
      fit_unit = "Selected weighted Step 2 case"))
  if (!is.null(regression)) model$regression <- regression
  model$marginals <- query_bn_static_cases(list(model = model))$marginals
  model
}

query_prepared_bn_v1 <- function(request) {
  model <- request$model
  if (!identical(model$prepared_schema_version, "prepared_bn_model_v1"))
    pa_fail("Expected a frozen prepared BN model")
  graph <- bs_validate_model(model)
  evidence <- pa_value(request$evidence, list())
  if (length(evidence) && (is.null(names(evidence)) ||
      any(!names(evidence) %in% graph$ids))) pa_fail("Evidence must name known nodes")
  answer <- query_bn_static_cases(list(model = model, evidence = evidence,
    targets = pa_value(request$targets, as.list(graph$ids))))
  built <- bs_grain(model, graph)
  if (length(evidence)) {
    engine <- gRain::setEvidence(built$engine, nodes = unname(built$safe[names(evidence)]),
      states = bs_vector(evidence))
    likelihood <- gRain::pEvidence(engine)
  } else likelihood <- 1
  if (!is.finite(likelihood) || likelihood <= 0)
    pa_fail("Impossible evidence in the selected frozen model")
  answer$evidence_probability <- as.numeric(likelihood)
  answer$view_digest <- model$view_digest
  answer
}

pa_collection_signature <- function(model) {
  graph <- bs_graph(model$bn)
  ids <- sort(graph$ids)
  nodes <- setNames(lapply(ids, function(id) {
    node <- graph$nodes[[match(id, graph$ids)]]
    time <- pa_value(node$fit_time, NULL)
    list(node_id = as.character(pa_value(node$node_id, id)),
      fit_time = if (is.null(time)) NULL else as.numeric(unlist(time, use.names = FALSE)))
  }), ids)
  edge_value <- function(edge, key) {
    value <- edge[[key]]
    if (is.null(value)) "" else paste(unlist(value, use.names = FALSE), collapse = ":")
  }
  edges <- sort(vapply(graph$edges, function(edge)
    paste(edge_value(edge, "from"), edge_value(edge, "to"),
      edge_value(edge, "kind"), edge_value(edge, "delta_time"), sep = "|"), character(1)))
  states <- setNames(lapply(ids, function(id) bs_vector(model$bn$node_states[[id]])), ids)
  information <- model$discretization$nodes
  information_ids <- vapply(information, function(node)
    as.character(node$node_id), character(1))
  if (length(information_ids) != length(ids) || anyDuplicated(information_ids) ||
      !setequal(information_ids, ids))
    pa_fail("BN collection requires one state definition for every graph node")
  cuts <- setNames(lapply(ids, function(id) {
    info <- information[[match(id, information_ids)]]
    list(base_node_id = as.character(pa_value(info$base_node_id, id)),
      n_states = as.integer(unlist(info$n_states, use.names = FALSE)),
      cutpoints = as.numeric(unlist(info$cutpoints, use.names = FALSE)))
  }), ids)
  list(nodes = nodes, edges = edges, states = states, cuts = cuts)
}

query_prepared_bn_collection_v1 <- function(request) {
  models <- request$models
  if (!is.list(models) || !length(models) || length(models) > 100L)
    pa_fail("A bounded nonempty collection of complete BN models is required")
  reference <- pa_collection_signature(models[[1]])
  states <- models[[1]]$bn$node_states
  graph_ids <- bs_graph(models[[1]]$bn)$ids
  for (model in models) {
    if (!identical(model$prepared_schema_version, "prepared_bn_model_v1") ||
        !identical(model$view_digest, models[[1]]$view_digest) ||
        !identical(pa_collection_signature(model), reference))
      pa_fail("BN collection components require the same graph and frozen state basis")
  }
  weights <- pa_value(request$component_weights, rep(1 / length(models), length(models)))
  if (!is.numeric(weights) || length(weights) != length(models) ||
      any(!is.finite(weights)) || any(weights < 0) || sum(weights) <= 0)
    pa_fail("Invalid complete-model component weights")
  weights <- weights / sum(weights)
  evidence <- pa_value(request$evidence, list())
  targets <- pa_value(request$targets, as.list(graph_ids))
  answers <- lapply(models, function(model) tryCatch(
    query_prepared_bn_v1(list(model = model, evidence = evidence, targets = targets)),
    error = function(e) if (grepl("Impossible evidence", conditionMessage(e), fixed = TRUE)) NULL else stop(e)))
  likelihoods <- vapply(answers, function(x) if (is.null(x)) 0 else x$evidence_probability, numeric(1))
  fixed_exposure <- "Exposure" %in% graph_ids && "Exposure" %in% names(evidence)
  if (fixed_exposure) {
    exposure_condition <- list(Exposure = evidence$Exposure)
    exposure_likelihoods <- vapply(models, function(model) {
      conditioned <- tryCatch(query_prepared_bn_v1(list(model = model,
        evidence = exposure_condition, targets = "Exposure")),
        error = function(e) if (grepl("Impossible evidence", conditionMessage(e), fixed = TRUE)) NULL else stop(e))
      if (is.null(conditioned)) 0 else conditioned$evidence_probability
    }, numeric(1))
    if (any(exposure_likelihoods <= 0 & weights > 0))
      pa_fail("Fixed Exposure state has no support in one or more collection components")
    likelihoods <- ifelse(exposure_likelihoods > 0,
      likelihoods / pmax(exposure_likelihoods, .Machine$double.xmin), 0)
  }
  mass <- weights * likelihoods
  if (sum(mass) <= 0) pa_fail("Impossible evidence in every complete BN model")
  posterior <- mass / sum(mass)
  targets <- bs_vector(targets)
  marginals <- setNames(lapply(targets, function(target) {
    labels <- bs_vector(states[[target]])
    p <- vapply(labels, function(label) sum(vapply(seq_along(answers), function(i)
      if (is.null(answers[[i]])) 0 else posterior[i] * answers[[i]]$marginals[[target]][[label]], numeric(1))), numeric(1))
    setNames(as.list(p), labels)
  }), targets)
  list(status = "success", marginals = marginals,
    prior_component_weights = as.list(weights),
    posterior_component_weights = as.list(posterior),
    component_evidence_likelihoods = as.list(likelihoods),
    interpretation = if (fixed_exposure)
      "Mixture weights update from biological evidence conditional on fixed Exposure; CPTs were not averaged" else
      "Evidence-reweighted mixture of complete frozen networks; CPTs were not averaged")
}
