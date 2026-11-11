# Refit the three archived tPSEM and DBN demo pairing realizations.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(paste(
    "Usage: Rscript --vanilla r/prepared/run.R",
    "PREPARED_CASES_CSV TEMPORAL_EDGES_CSV FOLD_CASES_CSV OUTPUT_ROOT"
  ), call. = FALSE)
}

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
source(file.path(script_dir, "bn_query_core.R"), local = TRUE)
source(file.path(script_dir, "prepared_fit_core.R"), local = TRUE)

RNGkind("Mersenne-Twister", "Inversion", "Rejection")

read_cases <- function(path) {
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
    fileEncoding = "UTF-8-BOM")
}

write_result <- function(path, value) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(value, path, auto_unbox = TRUE, digits = NA,
    pretty = TRUE, null = "null", na = "null")
}

write_details <- function(directory, model, cases) {
  write_result(file.path(directory, "frozen_model.json"), model)
  utils::write.csv(cases, file.path(directory, "prepared_cases.csv"),
    row.names = FALSE, na = "")
}

full_cases <- read_cases(args[1L])
fold_cases <- read_cases(args[3L])
out_dir <- args[4L]
base_ids <- c("KE1446", "KE1521", "KE1771", "KE1821")
time_points <- c(2, 7, 2, 4)
node_ids <- paste0(base_ids, "@", time_points)
node_labels <- c("Decrease, Coupling of OXPHOS", "Decrease, Growth",
  "Decrease, ATP pool", "Decrease, Cell proliferation")
nodes <- lapply(seq_along(node_ids), function(i) {
  list(id = node_ids[[i]], node_id = base_ids[[i]], fit_time = time_points[[i]],
    label = node_labels[[i]], role = "KE")
})
edge_table <- read_cases(args[2L])
expected_edges <- c("KE1446@2->KE1771@2", "KE1771@2->KE1821@4",
  "KE1821@4->KE1521@7")
if (!identical(paste(edge_table$from, edge_table$to, sep = "->"), expected_edges) ||
    !all(edge_table$kind == "aop") ||
    !isTRUE(all.equal(edge_table$delta_time, c(0, 2, 3)))) {
  stop("The prepared edges do not match the demo path.", call. = FALSE)
}
network <- list(nodes = nodes, edges = bs_records(edge_table))
case_md5 <- unname(tools::md5sum(args[1L]))
manifest <- list(schema_version = "step2_prepared_view_v1",
  generation = "lemna-cccp-v16-local", view_digest = case_md5,
  view_kind = "constructed_cases")

run_ids <- as.character(0:2)
realizations <- paste0(run_ids, ":0")
if (!identical(sort(unique(full_cases$realization_id)), realizations) ||
    !all(full_cases$time_context_id == "time_context_1") ||
    !all(table(full_cases$realization_id) == 30L)) {
  stop("Expected three 30-case demo realizations.", call. = FALSE)
}

fit_request <- function(model_type, cases, config = list()) {
  list(schema_version = "prepared_analysis_request_v1", model_type = model_type,
    prepared_view = list(manifest = manifest, cases = cases),
    network = network, config = config)
}

context_identity <- function(index) {
  list(time_context_id = "time_context_1", run_id = run_ids[[index]],
    realization_id = realizations[[index]])
}

context_result <- function(index, cases, model) {
  result <- list(context = context_identity(index), case_count = nrow(cases),
    weight_sum = sum(cases$case_weight), network = network)
  if (identical(model$schema_version, "prepared_psem_model_v1")) {
    result$equations <- model$equations
    result$network$edges <- lapply(network$edges, function(edge) {
      equation <- model$equations[[edge$to]]
      edge$coefficient <- equation$coefficients[[edge$from]]
      edge$coefficient_kind <- "raw_weighted_component"
      edge$path_inference <- equation$path_inference[[edge$from]]
      edge
    })
  } else {
    result$bn <- model$bn
    result$discretization <- model$discretization
    result$marginals <- model$marginals
  }
  result
}

context_cases <- function(cases, index) {
  cases[cases$realization_id == realizations[[index]], , drop = FALSE]
}

run_directory <- function(directory, index) {
  file.path(directory, paste0("realization-", run_ids[[index]], "-0"))
}

assess_psem <- function(cases, index) {
  target <- "KE1521@7"
  treatments <- sort(unique(cases$condition_id[cases$fit_dose != 0]))
  partition_cases <- context_cases(fold_cases, index)
  partition_cases <- partition_cases[partition_cases$held_out_condition_id %in% treatments, ]
  if (length(treatments) != 5L ||
      !setequal(treatments, partition_cases$held_out_condition_id) ||
      !setequal(c("training", "test"), partition_cases$partition)) {
    stop("Expected one source-level fold for each of five treated doses.", call. = FALSE)
  }
  folds <- lapply(treatments, function(condition) {
    fold <- partition_cases[partition_cases$held_out_condition_id == condition, ]
    training <- fold[fold$partition == "training", names(cases)]
    test <- fold[fold$partition == "test", names(cases)]
    if (nrow(training) != 25L || nrow(test) != 5L ||
        any(training$condition_id == condition) || any(training$fit_dose < 0) ||
        !any(training$fit_dose == 0) || !all(test$condition_id == condition)) {
      stop("Fold partitions do not match the demo design.", call. = FALSE)
    }
    model <- run_prepared_psem_v1(fit_request("psem", training))
    score <- score_prepared_psem_v1(list(model = model, cases = test,
      target = target, mode = "one_step"))
    condition_means <- vapply(split(training, training$condition_id),
      function(rows) weighted.mean(rows[[target]], rows$case_weight), numeric(1L))
    reference_residuals <- test[[target]] - mean(condition_means)
    reference <- list(
      rmse = sqrt(weighted.mean(reference_residuals^2, test$case_weight)),
      mae = weighted.mean(abs(reference_residuals), test$case_weight)
    )
    list(held_out_condition_id = condition, training_cases = nrow(training),
      test_cases = nrow(test), score = score, reference = reference)
  })
  metric <- function(section, name) {
    vapply(folds, function(fold) fold[[section]][[name]], numeric(1L))
  }
  list(status = "success", target = target, realization_id = realizations[[index]],
    split_unit = "held_out_treatment_condition", folds = folds, summary = list(
      fold_count = length(folds), scored_cases = sum(vapply(folds, function(fold)
        fold$test_cases, integer(1L))),
      rmse = sqrt(mean(metric("score", "rmse")^2)), mae = mean(metric("score", "mae")),
      reference_model = "training_mean_response",
      reference_rmse = sqrt(mean(metric("reference", "rmse")^2)),
      reference_mae = mean(metric("reference", "mae"))))
}

run_psem <- function() {
  directory <- file.path(out_dir, "tPSEM")
  contexts <- lapply(seq_along(realizations), function(index) {
    cases <- context_cases(full_cases, index)
    model <- run_prepared_psem_v1(fit_request("psem", cases))
    destination <- run_directory(directory, index)
    write_details(destination, model, cases)
    write_result(file.path(destination, "assessment", "results.json"),
      assess_psem(cases, index))
    context_result(index, cases, model)
  })
  write_result(file.path(directory, "results.json"),
    list(status = "success", contexts = contexts))
}

run_dbn <- function(destination, mode, method, state_count) {
  directory <- file.path(out_dir, destination)
  basis <- NULL
  models <- lapply(seq_along(realizations), function(index) {
    cases <- context_cases(full_cases, index)
    design <- setNames(lapply(base_ids, function(id)
      list(method = method, n_states = state_count)), base_ids)
    if (!is.null(basis)) {
      for (state in basis) {
        design[[state$base_node_id]]$fixed_cutpoints <- state$cutpoints
      }
    }
    config <- list(cpt_mode = mode, smoothing_alpha = 1,
      state_design = design, regression = list(simulation = list(n = 10000, seed = 42)))
    model <- run_prepared_bn_v1(fit_request("bn", cases, config))
    if (is.null(basis)) basis <<- model$discretization$nodes
    write_details(run_directory(directory, index), model, cases)
    model
  })
  contexts <- lapply(seq_along(models), function(index)
    context_result(index, context_cases(full_cases, index), models[[index]]))
  write_result(file.path(directory, "results.json"),
    list(status = "success", contexts = contexts))

  evidence <- list(`KE1446@2` = "state_1")
  targets <- as.list(c("KE1446@2", "KE1771@2", "KE1821@4", "KE1521@7"))
  write_result(file.path(directory, "baseline", "results.json"),
    query_prepared_bn_collection_v1(list(models = models, targets = targets)))
  write_result(file.path(directory, "query.json"),
    query_prepared_bn_collection_v1(list(models = models, evidence = evidence,
      targets = targets)))
  for (index in seq_along(models)) {
    write_result(file.path(run_directory(directory, index), "query", "results.json"),
      query_prepared_bn_v1(list(model = models[[index]], evidence = evidence,
        targets = targets)))
  }
}

run_psem()
run_dbn("DBN", "bayesian_regression", "quantile", 3L)
run_dbn(file.path("DBN", "comparison", "counts-3"), "observed_counts", "quantile", 3L)
run_dbn(file.path("DBN", "comparison", "residual-3"), "regression", "quantile", 3L)
run_dbn(file.path("DBN", "comparison", "residual-5"), "regression", "equal_width", 5L)
