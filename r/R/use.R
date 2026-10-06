# Taking the numbers out of the object, where the caveats are acted on rather than explained.

#' The consumed and produced matrices
#'
#' Taxa as rows, resources (metabolites) as columns, each cell the amount the taxon removed or produced in
#' its monoculture over the phase, in mM (or 1 and 0 when foodnet reported booleans). `NA` is never zero:
#' it is a compound that was not assayed for that taxon, or a link seen only in another medium. The
#' evidence matrices (`x$evidence_consumed`, `x$evidence_produced`) say which.
#'
#' @param x CRM parameters from [foodnet_listen()] or [foodnet_crm()].
#' @return A numeric matrix with taxa and resources as its row and column names.
#' @export
crm_consumed <- function(x) {
    stopifnot(inherits(x, "foodnet_crm"))
    x$consumed
}

#' @rdname crm_consumed
#' @export
crm_produced <- function(x) {
    stopifnot(inherits(x, "foodnet_crm"))
    x$produced
}

#' The growth rates
#'
#' Each taxon's maximum specific growth rate in monoculture, as foodnet derived it: from the replicates
#' whose metabolites gave the values, or else from another monoculture in the same medium. A taxon with
#' none is `NA`.
#'
#' @param x CRM parameters from [foodnet_listen()] or [foodnet_crm()].
#' @param missing A number to use for the taxa that have no rate, or `NULL` to leave them `NA`. Giving
#'   one is a choice about data that is not there, so it is never made for you.
#' @return A named numeric vector in the row order of the matrices.
#' @export
crm_rates <- function(x, missing = NULL) {
    stopifnot(inherits(x, "foodnet_crm"))
    rates <- x$growth_rates
    if (!is.null(missing)) {
        if (!is.numeric(missing) || length(missing) != 1) stop_foodnet("missing must be one number")
        rates[is.na(rates)] <- missing
    }
    rates
}

#' The initial resource concentrations
#'
#' Each metabolite's concentration at the first sample of the cultures the values come from, averaged,
#' in mM: the medium a simulation starts from. A resource never assayed in that medium is `NA`.
#'
#' @param x CRM parameters from [foodnet_listen()] or [foodnet_crm()].
#' @param missing A concentration for the resources that have none, or `NULL` to leave them `NA`.
#' @return A named numeric vector in the column order of the matrices.
#' @export
crm_resources <- function(x, missing = NULL) {
    stopifnot(inherits(x, "foodnet_crm"))
    initial <- x$initial
    if (!is.null(missing)) initial[is.na(initial)] <- missing
    initial
}

#' The efficiency matrix of a consumer-resource model
#'
#' A consumer-resource model such as miaSim's `simulateConsumerResource` takes an efficiency matrix `E`,
#' taxa by resources, positive where a taxon consumes a resource and negative where it produces one. This
#' builds one from the measured amounts.
#'
#' With `scale = "miasim"` (the default), `E` follows miaSim's own equations, read from its source
#' (`consumerResourceModel`, miaSim 1.18): a taxon grows by `growth_rate * sum_j E[j] * R[j] / (R[j] + K[j])`
#' per unit of its abundance, takes up each resource it consumes at `R / (R + K)` per unit of abundance
#' whatever the size of `E`, and makes each by-product at `|E|` times its growth. So a positive entry is a
#' yield, the biomass made per mM taken up, and a negative one the mM of a by-product made per unit of
#' biomass divided by the growth rate. For each taxon, with its biomass change `dx` and growth rate `mu`
#' over the phase, `E` is `dx / (mu * total consumed)` for every resource it consumed and
#' `-produced * mu / dx` for every resource it produced, so that a taxon simulated alone gains its measured
#' biomass from its measured uptake and makes its measured by-products. How its uptake splits between
#' resources is set by the Monod constants, which foodnet does not measure; [crm_backcheck()] says how close
#' a choice of them comes. The biomass is in the unit of each taxon's growth curve (`x$biomass_unit`), and so
#' must be its starting abundance in a simulation.
#'
#' `scale = "shares"` is foodnet 0.1.0's matrix: each row divided by the total consumed, so consumed entries
#' are shares of uptake and produced entries by-product per unit taken up. miaSim does not read `E` that way
#' (it reads a positive entry as a yield), so a simulation from it does not reproduce the monocultures; it is
#' kept for other models. `scale = "none"` is consumed minus produced, in mM.
#'
#' `NA` cells have no size. With `na = "zero"` they count as 0, and a warning says how many of each kind
#' were set to 0: never assayed, seen only in another medium (a real link), inconclusive (measured, but its
#' spread reaches across the limit) and without a phase. `"stop"` refuses until you fill them yourself.
#'
#' @param x CRM parameters from [foodnet_listen()] or [foodnet_crm()].
#' @param scale `"miasim"`, `"shares"` or `"none"`, as above.
#' @param na `"zero"` or `"stop"`, as above.
#' @return A numeric matrix, taxa by resources.
#' @examples
#' \dontrun{
#' E <- crm_efficiency(crm)
#' }
#' @export
crm_efficiency <- function(x, scale = c("miasim", "shares", "none"), na = c("zero", "stop")) {
    stopifnot(inherits(x, "foodnet_crm"))
    scale <- match.arg(scale)
    na <- match.arg(na)
    if (identical(x$values, "booleans") && scale != "none") {
        stop_foodnet("these parameters are booleans (1 = it happened), not amounts: an efficiency matrix ",
                     "needs amounts. Run the search again without Report everything as booleans.")
    }
    consumed <- x$consumed
    produced <- x$produced
    if (anyNA(consumed) || anyNA(produced)) {
        if (na == "stop") {
            stop_foodnet(sum(is.na(consumed)) + sum(is.na(produced)), " cell(s) are NA (not assayed, seen ",
                         "only in another medium, inconclusive or without a phase). Fill them, or use ",
                         "crm_efficiency(x, na = \"zero\").")
        }
        zeroed <- c(x$evidence_consumed[is.na(consumed)], x$evidence_produced[is.na(produced)])
        kinds <- table(zeroed)
        warning("NA cells counted as 0: ", paste0(kinds, " ", names(kinds), collapse = ", "),
                if ("presence_only" %in% names(kinds)) " (presence_only links are real; their size is unknown)",
                call. = FALSE)
        consumed[is.na(consumed)] <- 0
        produced[is.na(produced)] <- 0
    }
    elsewhere <- sum(x$evidence_consumed == "presence_only" & consumed > 0, na.rm = TRUE) +
        sum(x$evidence_produced == "presence_only" & produced > 0, na.rm = TRUE)
    if (elsewhere) {
        warning(elsewhere, " amount(s) were measured in another medium than the values (foodnet's \"value\" ",
                "setting for links seen only there) and enter E as they are", call. = FALSE)
    }
    if (scale == "none") return(consumed - produced)
    total <- rowSums(consumed)
    none <- names(total)[total <= 0]
    if (length(none)) {
        warning("no measured uptake for ", paste(none, collapse = ", "), ": their consumption rows are 0, so a ",
                "CRM gives them nothing to grow on", call. = FALSE)
    }
    if (scale == "shares") {
        return((consumed - produced) / ifelse(total > 0, total, 1))
    }
    rates <- x$growth_rates
    dx <- x$biomass_change
    lacking <- names(rates)[is.na(rates) | is.na(dx) | dx <= 0]
    if (length(lacking)) {
        stop_foodnet("miaSim's yields need a growth rate and a positive biomass change for every taxon; ",
                     paste(lacking, collapse = ", "), " lack one (see print(x)). Leave them out, or use ",
                     "crm_efficiency(x, scale = \"shares\") for a model that reads E as shares.")
    }
    yield <- ifelse(total > 0, dx / (rates * total), 0)
    (consumed > 0) * yield - produced * (rates / dx)
}

#' The README foodnet wrote with these numbers
#'
#' @param x CRM parameters from [foodnet_listen()] or [foodnet_crm()].
#' @param print Set to `FALSE` to return the text without printing it.
#' @return The text, invisibly when printed.
#' @export
crm_readme <- function(x, print = TRUE) {
    stopifnot(inherits(x, "foodnet_crm"))
    if (print) {
        cat(x$readme)
        return(invisible(x$readme))
    }
    x$readme
}

#' Write the parameters as files
#'
#' `consumed.csv`, `produced.csv`, `growth_rates.csv`, `initial_concentrations.csv` and `README.txt`, as
#' foodnet's download holds them, so what arrived over the wire can be kept beside the analysis.
#'
#' @param x CRM parameters from [foodnet_listen()] or [foodnet_crm()].
#' @param dir Directory to write into. It is created when it does not exist.
#' @return The paths written, invisibly.
#' @export
crm_write <- function(x, dir) {
    stopifnot(inherits(x, "foodnet_crm"))
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
    paths <- file.path(dir, c("consumed.csv", "produced.csv", "growth_rates.csv",
                              "initial_concentrations.csv", "README.txt"))
    write.csv(x$consumed, paths[1])
    write.csv(x$produced, paths[2])
    write.csv(data.frame(taxon = x$taxa, growth_rate = unname(x$growth_rates), unit = x$growth_rate_unit,
                         stringsAsFactors = FALSE), paths[3], row.names = FALSE)
    write.csv(data.frame(resource = x$resources, initial_mM = unname(x$initial), stringsAsFactors = FALSE),
              paths[4], row.names = FALSE)
    writeLines(x$readme, paths[5])
    invisible(paths)
}

#' Shape the parameters for miaSim
#'
#' Returns the arguments `miaSim::simulateConsumerResource` takes: the numbers of taxa and resources and
#' their names, the efficiency matrix `E`, the initial resource concentrations, the growth rates, the
#' starting abundances and the Monod constants. miaSim is not needed to call this, and no other simulator
#' is assumed: the result is a plain list, so
#' `do.call(miaSim::simulateConsumerResource, c(as_miasim(crm, x0 = ..., monod_constant = ...), list(t_end = 48, t_store = 480)))`
#' runs it. miaSim stores `t_store` of its steps of 0.1 h, evenly spread; with `t_store = 10 * t_end` it
#' keeps every step and the run reaches `t_end` (fewer, and it stops early).
#'
#' A simulation needs a growth rate for every taxon, a starting concentration for every resource, a
#' starting abundance for every taxon and Monod constants, so this stops when one is missing rather than
#' passing a number nobody measured: miaSim would otherwise draw starting abundances and Monod constants at
#' random. It also stops for parameters that pool every medium (their starting concentrations describe no
#' real medium) or come from the stationary phase (uptake without growth), unless `allow` names them.
#'
#' @param x CRM parameters from [foodnet_listen()] or [foodnet_crm()].
#' @param x0 Starting abundances of the taxa, one per taxon in the order of `x$taxa` (or named), each in the
#'   unit of that taxon's growth curve (`x$biomass_unit`); `x$biomass_start` is each taxon's own start.
#' @param monod_constant Monod constants in mM: one number for all, or a taxa by resources matrix.
#' @param E An efficiency matrix, by default [crm_efficiency()] of `x`.
#' @param missing_rate A growth rate for the taxa that have none, or `NULL` to stop when any does.
#' @param missing_resource A starting concentration for the resources that have none, or `NULL` to stop.
#' @param allow Caveats to accept: `"mixed_media"`, `"stationary_phase"`.
#' @return A list with `n_species`, `n_resources`, `names_species`, `names_resources`, `E`, `x0`,
#'   `resources`, `growth_rates` and `monod_constant`.
#' @examples
#' \dontrun{
#' args <- as_miasim(crm, x0 = crm$biomass_start, monod_constant = 1)
#' set.seed(1)
#' tse <- do.call(miaSim::simulateConsumerResource, c(args, list(t_end = 48, t_store = 480)))
#' }
#' @export
as_miasim <- function(x, x0, monod_constant, E = NULL, missing_rate = NULL, missing_resource = NULL,
                      allow = character()) {
    stopifnot(inherits(x, "foodnet_crm"))
    if (x$caveats$mixed_media && !"mixed_media" %in% allow) {
        stop_foodnet("these values pool every medium (Ignore media differences): their starting concentrations ",
                     "describe no real medium. Run the search with one medium, or pass allow = \"mixed_media\".")
    }
    if (x$caveats$stationary_phase && !"stationary_phase" %in% allow) {
        stop_foodnet("these are stationary-phase amounts: uptake without growth, which a CRM reads as growth. ",
                     "Use the exponential phase, or pass allow = \"stationary_phase\".")
    }
    if (missing(x0) || is.null(x0)) {
        stop_foodnet("give x0, a starting abundance per taxon in the unit of its growth curve (x$biomass_unit); ",
                     "x$biomass_start holds each taxon's own start. miaSim would otherwise draw them at random.")
    }
    if (missing(monod_constant) || is.null(monod_constant)) {
        stop_foodnet("give monod_constant, in mM (one number, or a taxa by resources matrix); foodnet measures ",
                     "none, and miaSim would otherwise draw them at random. crm_backcheck() tests a choice.")
    }
    n <- length(x$taxa)
    m <- length(x$resources)
    x0 <- if (!is.null(names(x0)) && all(x$taxa %in% names(x0))) x0[x$taxa] else x0
    if (length(x0) != n || anyNA(x0)) stop_foodnet("x0 needs one number per taxon (", n, ")")
    K <- if (length(monod_constant) == 1) matrix(monod_constant, n, m) else as.matrix(monod_constant)
    if (!identical(dim(K), c(n, m)) || anyNA(K)) stop_foodnet("monod_constant needs one number or a ", n, " by ",
                                                             m, " matrix")
    rates <- crm_rates(x, missing = missing_rate)
    if (anyNA(rates)) {
        stop_foodnet("no growth rate for ", paste(names(rates)[is.na(rates)], collapse = ", "),
                     ": a CRM simulation needs one for every taxon. Give as_miasim(x, missing_rate = ) a ",
                     "number of your own, or leave those taxa out.")
    }
    resources <- crm_resources(x, missing = missing_resource)
    if (anyNA(resources)) {
        stop_foodnet("no starting concentration for ", paste(names(resources)[is.na(resources)], collapse = ", "),
                     ": give as_miasim(x, missing_resource = ) a number (0 for a resource the medium lacks).")
    }
    if (is.null(E)) E <- crm_efficiency(x)
    list(n_species = n, n_resources = m, names_species = x$taxa, names_resources = x$resources, E = E,
         x0 = unname(as.numeric(x0)), resources = unname(resources), growth_rates = unname(rates),
         monod_constant = unname(K))
}

#' Simulate each taxon alone and compare with its monoculture
#'
#' Runs miaSim's consumer-resource model for each taxon on its own, from its measured starting abundance
#' and the medium's starting concentrations, over the length of its phase, and compares what it gains and
#' what it takes up and makes with what its monoculture did. With [crm_efficiency()]'s default scale the
#' biomass and the by-products follow the measurements by construction whenever the total uptake does; the
#' uptake of each resource follows the Monod constants, so this is where a choice of them is tested.
#'
#' @param x CRM parameters from [foodnet_listen()] or [foodnet_crm()].
#' @param monod_constant Monod constants in mM: one number, or a taxa by resources matrix.
#' @param E An efficiency matrix, by default [crm_efficiency()] of `x`.
#' @param x0 Starting abundances, by default each taxon's own (`x$biomass_start`).
#' @param missing_resource A starting concentration for resources with none (default 0).
#' @return A data frame: taxon, what (`"biomass"`, or a resource), direction, measured, simulated and their
#'   ratio; printed in short.
#' @examples
#' \dontrun{
#' crm_backcheck(crm, monod_constant = 1)
#' }
#' @export
crm_backcheck <- function(x, monod_constant, E = NULL, x0 = x$biomass_start, missing_resource = 0) {
    stopifnot(inherits(x, "foodnet_crm"))
    if (!requireNamespace("miaSim", quietly = TRUE)) {
        stop_foodnet("crm_backcheck() runs miaSim: install it with BiocManager::install(\"miaSim\")")
    }
    if (missing(monod_constant)) stop_foodnet("give monod_constant, in mM; foodnet measures none")
    if (is.null(E)) E <- crm_efficiency(x)
    n <- length(x$taxa)
    m <- length(x$resources)
    K <- if (length(monod_constant) == 1) matrix(monod_constant, n, m) else as.matrix(monod_constant)
    start <- crm_resources(x, missing = missing_resource)
    consumed <- x$consumed
    produced <- x$produced
    rows <- list()
    for (i in seq_len(n)) {
        taxon <- x$taxa[i]
        hours <- x$phase_hours[[taxon]]
        if (is.na(hours) || is.na(x0[i]) || is.na(x$growth_rates[[taxon]])) next
        steps <- max(1L, round(hours / 0.1))
        tse <- miaSim::simulateConsumerResource(
            n_species = 1, n_resources = m, names_species = taxon, names_resources = x$resources,
            E = E[i, , drop = FALSE], x0 = x0[[i]], resources = unname(start), growth_rates = x$growth_rates[[taxon]],
            monod_constant = K[i, , drop = FALSE], t_end = (steps + 1) * 0.1, t_store = steps + 1)
        abundance <- SummarizedExperiment::assay(tse)[1, ]
        times <- SummarizedExperiment::colData(tse)$time
        last <- which.min(abs(times - hours))
        left <- S4Vectors::metadata(tse)$resources[last, seq_len(m)]
        change <- left - start
        rows[[length(rows) + 1]] <- data.frame(taxon = taxon, what = "biomass", direction = "grown",
                                               measured = x$biomass_change[[taxon]],
                                               simulated = abundance[last] - abundance[1])
        for (j in seq_len(m)) {
            for (direction in c("consumed", "produced")) {
                measured <- if (direction == "consumed") consumed[i, j] else produced[i, j]
                simulated <- if (direction == "consumed") max(0, -change[j]) else max(0, change[j])
                if ((is.na(measured) || measured == 0) && simulated < 1e-6) next
                rows[[length(rows) + 1]] <- data.frame(taxon = taxon, what = x$resources[j], direction = direction,
                                                       measured = measured, simulated = simulated)
            }
        }
    }
    out <- do.call(rbind, rows)
    if (is.null(out)) stop_foodnet("no taxon has a growth rate, a starting abundance and a phase length")
    out$ratio <- out$simulated / out$measured
    rownames(out) <- NULL
    out
}

