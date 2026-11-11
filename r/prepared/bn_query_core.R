# Functions copied from qAOPweaver/qAOPweaverR/R/bn_static_cases.R
# qAOPweaver v0.29.0 revision 0143dd56dd21a5268b17ecbe8383647b1c6323c9.
# See the parent repository's provenance and licensing documentation.

bs_default <- function(x, default) if (is.null(x)) default else x

bs_fail <- function(...) stop(paste0(...), call. = FALSE)

bs_vector <- function(x) as.character(unlist(x, use.names = FALSE))

bs_records <- function(x) {
  if (is.data.frame(x)) return(lapply(seq_len(nrow(x)), function(i) as.list(x[i, , drop = FALSE])))
  bs_default(x, list())
}

bs_frame <- function(x) {
  if (is.data.frame(x)) return(x)
  if (!length(x)) return(data.frame())
  # Keep absent/null values as NA, not omitted columns.
  fields <- unique(unlist(lapply(x, names), use.names = FALSE))
  rows <- lapply(x, function(r) setNames(lapply(fields, function(f) bs_default(r[[f]], NA)), fields))
  do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE, check.names = FALSE))
}

bs_number <- function(x, default, name, lower = -Inf, upper = Inf, integer = FALSE) {
  x <- bs_default(x, default)
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < lower || x > upper ||
      (integer && x != floor(x))) bs_fail(name, " must be finite and between ", lower, " and ", upper)
  as.numeric(x)
}

bs_graph <- function(network) {
  nodes <- bs_records(network$nodes); edges <- bs_records(network$edges)
  ids <- vapply(nodes, function(n) as.character(bs_default(n$id, "")), character(1))
  if (!length(ids) || length(ids) > 30 || any(!nzchar(ids)) || anyDuplicated(ids))
    bs_fail("Graph requires 1-30 unique nonempty node IDs")
  parents <- setNames(rep(list(character()), length(ids)), ids)
  for (e in edges) {
    if (length(e$from) != 1L || length(e$to) != 1L || !e$from %in% ids || !e$to %in% ids)
      bs_fail("Every edge must reference declared nodes")
    if (e$from == e$to || e$from %in% parents[[e$to]]) bs_fail("Self edges or duplicate edges are not allowed")
    parents[[e$to]] <- c(parents[[e$to]], e$from)
  }
  remaining <- ids; order <- character()
  while (length(remaining)) {
    ready <- remaining[vapply(remaining, function(n) all(parents[[n]] %in% order), logical(1))]
    if (!length(ready)) bs_fail("Graph must be acyclic: cycle detected")
    order <- c(order, ready); remaining <- setdiff(remaining, ready)
  }
  list(nodes = nodes, edges = edges, ids = ids, parents = parents, order = order)
}

bs_limits <- function(graph, states) {
  sizes <- vapply(states, length, integer(1))
  cells <- vapply(graph$ids, function(n) prod(sizes[c(n, graph$parents[[n]])]), numeric(1))
  if (any(cells > 50000) || sum(cells) > 250000) bs_fail("CPT cell budget exceeded (50000/node, 250000 total)")
  # Bound induced factor size before exact inference, without joint enumeration.
  adjacent <- setNames(rep(list(character()), length(graph$ids)), graph$ids)
  for (n in graph$ids) {
    family <- c(n, graph$parents[[n]])
    for (v in family) adjacent[[v]] <- union(adjacent[[v]], setdiff(family, v))
  }
  active <- graph$ids
  while (length(active)) {
    factor_sizes <- vapply(active, function(n) prod(sizes[c(n, intersect(adjacent[[n]], active))]), numeric(1))
    n <- active[which.min(factor_sizes)]
    if (min(factor_sizes) > 1000000) bs_fail("Exact-inference factor budget exceeded (1000000 cells)")
    neighbours <- intersect(adjacent[[n]], active)
    for (v in neighbours) adjacent[[v]] <- union(adjacent[[v]], setdiff(neighbours, v))
    active <- setdiff(active, n)
  }
  invisible(list(cells = cells, triangulated_adjacency = adjacent))
}

bs_assign <- function(x, cuts) paste0("state_", findInterval(x, cuts) + 1L)

bs_assignments <- function(parents, states) {
  if (!length(parents)) return(list(structure(list(), names = character())))
  grid <- do.call(expand.grid, c(lapply(states[parents], bs_vector), stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE))
  bs_records(grid)
}

bs_validate_model <- function(model) {
  if (!identical(model$schema_version, "bn_static_cases_model_v1")) bs_fail("Unsupported frozen model schema")
  graph <- bs_graph(model$bn); states <- model$bn$node_states
  if (!setequal(names(states), graph$ids) || !setequal(names(model$bn$cpts), graph$ids)) bs_fail("Frozen model requires states and one CPT for every node")
  for (n in graph$ids) {
    labels <- bs_vector(states[[n]])
    if (length(labels) < 2 || length(labels) > 10 || anyDuplicated(labels) || any(!nzchar(labels))) bs_fail("Invalid frozen state definitions for ", n)
  }
  bs_limits(graph, states)
  for (n in graph$ids) {
    cpt <- model$bn$cpts[[n]]; parents <- graph$parents[[n]]; labels <- bs_vector(states[[n]])
    if (!identical(bs_vector(cpt$parent_ids), parents)) bs_fail("CPT parent order differs from graph for ", n)
    expected <- bs_assignments(parents, states)
    if (length(cpt$rows) != length(expected)) bs_fail("Incomplete CPT configurations for ", n)
    seen <- character()
    for (r in cpt$rows) {
      if (!setequal(names(r$parent_assignment), parents)) bs_fail("Invalid CPT parent assignment for ", n)
      for (p in parents) if (!r$parent_assignment[[p]] %in% bs_vector(states[[p]])) bs_fail("Invalid parent state in CPT for ", n)
      key <- paste(vapply(parents, function(p) match(r$parent_assignment[[p]], bs_vector(states[[p]])), integer(1)), collapse = ":")
      if (key %in% seen) bs_fail("Duplicate CPT parent configuration for ", n)
      seen <- c(seen, key)
      if (!setequal(names(r$probabilities), labels)) bs_fail("CPT must name every child state for ", n)
      probabilities <- unlist(r$probabilities[labels], use.names = FALSE)
      if (!is.numeric(probabilities) || any(!is.finite(probabilities)) || any(probabilities < 0) || abs(sum(probabilities)-1) > 1e-8)
        bs_fail("CPT probabilities must be finite, nonnegative and sum to one for ", n)
    }
  }
  graph
}

bs_grain <- function(model, graph) {
  if (!requireNamespace("gRain", quietly = TRUE)) bs_fail("Exact static BN inference requires dependency gRain")
  safe <- setNames(paste0("N", seq_along(graph$ids)), graph$ids)
  cptables <- lapply(graph$ids, function(n) {
    parents <- graph$parents[[n]]; family <- c(n, parents)
    dims <- setNames(lapply(model$bn$node_states[family], bs_vector), safe[family])
    a <- array(NA_real_, dim = lengths(dims), dimnames = dims)
    for (r in model$bn$cpts[[n]]$rows) {
      index <- c(list(bs_vector(model$bn$node_states[[n]])), lapply(parents, function(p) r$parent_assignment[[p]]))
      a <- do.call(`[<-`, c(list(a), index, list(value = unlist(r$probabilities[bs_vector(model$bn$node_states[[n]])]))))
    }
    formula <- stats::as.formula(paste0("~", safe[[n]], if (length(parents)) paste0("|", paste(safe[parents], collapse = "+")) else ""))
    gRain::cptable(formula, values = as.vector(a), levels = bs_vector(model$bn$node_states[[n]]))
  })
  # Supply the same bounded triangulation to gRain, so its clique allocations
  # cannot choose an unrelated, larger elimination order after the guard.
  budget <- bs_limits(graph, model$bn$node_states)
  tug <- matrix(0, length(graph$ids), length(graph$ids), dimnames = list(unname(safe), unname(safe)))
  for (n in graph$ids) tug[safe[[n]], safe[budget$triangulated_adjacency[[n]]]] <- 1
  engine <- gRain::grain(gRain::compileCPT(cptables), compile = FALSE)
  list(engine = gRbase::compile(engine, tug = tug), safe = safe)
}

bs_query <- function(engine, safe, states, evidence, targets, allow_impossible = FALSE) {
  g <- engine
  if (length(evidence)) {
    g <- gRain::setEvidence(g, nodes = unname(safe[names(evidence)]), states = bs_vector(evidence))
    probability <- gRain::pEvidence(g)
    if (!is.finite(probability) || probability <= 0) {
      if (allow_impossible) return(NULL)
      bs_fail("Impossible evidence: probability is zero in the selected frozen model")
    }
  }
  raw <- gRain::querygrain(g, nodes = unname(safe[targets]), type = "marginal", exclude = FALSE)
  setNames(lapply(targets, function(n) {
    labels <- bs_vector(states[[n]]); p <- as.numeric(raw[[safe[[n]]]][labels])
    if (any(!is.finite(p)) || abs(sum(p)-1) > 1e-8) bs_fail("Exact inference returned invalid probabilities")
    setNames(as.list(p), labels)
  }), targets)
}

query_bn_static_cases <- function(request) {
  model <- request$model; graph <- bs_validate_model(model); evidence <- bs_default(request$evidence, list())
  if (length(evidence) && (is.null(names(evidence)) || anyDuplicated(names(evidence)) || any(!names(evidence) %in% graph$ids))) bs_fail("Evidence must name known nodes")
  for (n in names(evidence)) if (length(evidence[[n]]) != 1 || !evidence[[n]] %in% bs_vector(model$bn$node_states[[n]])) bs_fail("Unknown evidence state for ", n)
  targets <- bs_vector(bs_default(request$targets, as.list(graph$ids)))
  if (!length(targets) || any(!targets %in% graph$ids) || anyDuplicated(targets)) bs_fail("Targets must be unique known nodes")
  built <- bs_grain(model, graph)
  marginal <- bs_query(built$engine, built$safe, model$bn$node_states, evidence, targets)
  sensitivity <- list()
  if (isTRUE(request$sensitivity)) for (n in setdiff(graph$ids, names(evidence))) for (state in bs_vector(model$bn$node_states[[n]])) {
    changed <- evidence; changed[[n]] <- state
    posterior <- bs_query(built$engine, built$safe, model$bn$node_states, changed, targets, TRUE)
    sensitivity[[length(sensitivity)+1L]] <- list(evidence_node = n, evidence_state = state,
      status = if (is.null(posterior)) "impossible_evidence" else "success", marginals = posterior,
      maximum_absolute_change = if (is.null(posterior)) NULL else max(abs(unlist(posterior)-unlist(marginal))),
      interpretation = "Evidence-conditioned sensitivity on the same frozen CPTs; not an intervention")
  }
  list(status = "success", marginals = marginal, sensitivity = sensitivity, warnings = list())
}
