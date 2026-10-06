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
#' (`consumerResourceModel`, miaSim 1.18): a taxon of abundance `x` grows at
#' `growth_rate * x * sum_j E[j] * R[j] / (R[j] + K[j])`, takes up each resource it consumes at
#' `x * R / (R + K)` whatever the size of `E` (at most 1 mM per unit of abundance per hour), and makes each
#' by-product at `|E|` times its growth (without the growth rate). miaSim has no uptake rate of its own, so
#' its unit of abundance sets how fast a taxon eats: in cells/mL a culture would empty its medium within
#' minutes. foodnet therefore gives each taxon a unit of its own, [crm_scale()]: `dx * C / (mu * S)` units
#' of its growth curve, with `dx` its biomass change, `mu` its growth rate, `C` its total uptake and `S` the
#' sum of the squares of its uptakes. In that unit, `E` is each consumed resource's share of the uptake,
#' `c / C` (they sum to 1, so a saturated taxon grows at its measured rate, and a resource it took little of
#' weighs little), and `-produced * C / S` on each produced one, and a taxon alone gains its measured biomass
#' and makes its measured by-products in proportion to what it takes up, in about the measured time. The same data in another unit give the same simulation. How its uptake splits
#' between resources follows the Monod constants, which foodnet does not measure; [crm_backcheck()] says how
#' close a choice of them comes. [as_miasim()] puts the starting abundances into these units, and
#' [crm_unscale()] turns a simulation's abundances back into each growth curve's unit.
#'
#' `scale = "shares"` is foodnet 0.1.0's matrix: each row divided by the total consumed, so consumed entries
#' are shares of uptake and produced entries by-product per unit taken up; for other models.
#' `scale = "none"` is consumed minus produced, in mM.
#'
#' `NA` cells have no size, and foodnet never writes one as zero, so by default (`na = "stop"`) this refuses
#' until you fill them yourself or say how. With `na = "zero"` they count as 0, and a warning says how many
#' of each kind were set to 0: never assayed, seen only in another medium (a real link), inconclusive
#' (measured, but its replicates or experiments disagree) and without a phase.
#'
#' @param x CRM parameters from [foodnet_listen()] or [foodnet_crm()].
#' @param scale `"miasim"`, `"shares"` or `"none"`, as above.
#' @param na `"stop"` or `"zero"`, as above.
#' @param growth `"measured"` (the default): the growth rates as foodnet fitted them. `"phase_floor"`: at
#'   least the mean rate each taxon's own curves show over the phase (see [crm_scale()]).
#' @return A numeric matrix, taxa by resources.
#' @examples
#' \dontrun{
#' E <- crm_efficiency(crm, na = "zero")
#' }
#' @export
crm_efficiency <- function(x, scale = c("miasim", "shares", "none"), na = c("stop", "zero"),
                           growth = c("measured", "phase_floor")) {
    stopifnot(inherits(x, "foodnet_crm"))
    scale <- match.arg(scale)
    growth <- match.arg(growth)
    na <- match.arg(na)
    m <- amounts(x, na, booleans_ok = scale == "none")
    consumed <- m$consumed
    produced <- m$produced
    if (scale == "none") return(consumed - produced)
    total <- rowSums(consumed)
    none <- names(total)[total <= 0]
    if (length(none)) {
        warning("no measured uptake for ", paste(none, collapse = ", "), ": their consumption rows are 0, so a ",
                "CRM gives them nothing to grow on", if (scale == "miasim") ", and in miaSim nothing to make by-products from",
                call. = FALSE)
    }
    if (scale == "shares") {
        return((consumed - produced) / ifelse(total > 0, total, 1))
    }
    crm_scale(x, growth)        # stops, naming them, when a taxon lacks what the scale needs
    squares <- rowSums(consumed^2)
    consumed / ifelse(total > 0, total, 1) - produced * ifelse(squares > 0, total / squares, 0)
}

# The consumed and produced amounts with every NA set to 0 and counted by kind, or a stop.
#' @noRd
amounts <- function(x, na = "stop", booleans_ok = FALSE) {
    if (identical(x$values, "booleans") && !booleans_ok) {
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
    moving <- x$caveats$cautions
    moving <- if (is.null(moving)) character(0) else
        paste(moving$taxon, moving$resource)[grepl("still_changing", moving$cautions)]
    if (length(moving)) {
        # a script that never prints the parameters still hears it
        warning(length(moving), " value(s) miss use that went on after growth slowed (still_changing), so E ",
                "underrates them: ", paste(utils::head(moving, 4), collapse = ", "),
                if (length(moving) > 4) " and more" else "", ". A time window over both phases gives the whole change.",
                call. = FALSE)
    }
    elsewhere <- sum(x$evidence_consumed == "presence_only" & consumed > 0, na.rm = TRUE) +
        sum(x$evidence_produced == "presence_only" & produced > 0, na.rm = TRUE)
    if (elsewhere) {
        warning(elsewhere, " amount(s) were measured in another medium than the values (foodnet's \"value\" ",
                "setting for links seen only there) and enter E as they are", call. = FALSE)
    }
    list(consumed = consumed, produced = produced)
}

#' Each taxon's unit of abundance in a miaSim simulation
#'
#' How many units of a taxon's growth curve (`x$biomass_unit`) one unit of its abundance in a simulation
#' built by [crm_efficiency()] and [as_miasim()] stands for: `dx * C / (mu * S)`, with `dx` its biomass
#' change over the phase, `mu` its growth rate, `C` its total uptake in mM and `S` the sum of the squares of
#' its uptakes of each resource. In this unit miaSim's fixed uptake (1 mM per unit per hour at saturation) makes the taxon
#' grow at its measured rate with its measured yield. Divide starting abundances by it ([as_miasim()] does)
#' and multiply simulated ones by it ([crm_unscale()]).
#'
#' The growth rate in it is foodnet's fitted maximum rate, or with `growth = "phase_floor"` at least the
#' mean rate the taxon's own curves show over the phase. A maximum below that mean cannot be right: it comes
#' from a fitting window that took in the plateau of a fast, sparsely sampled curve (E. coli LF82: 0.375
#' fitted against 0.65 per hour over its phase), and it would grow the taxon too slowly. Such taxa are named in
#' a warning.
#'
#' @param x CRM parameters from [foodnet_listen()] or [foodnet_crm()].
#' @param growth `"measured"` or `"phase_floor"`, as above.
#' @return A named numeric vector, one per taxon.
#' @export
crm_scale <- function(x, growth = c("measured", "phase_floor")) {
    stopifnot(inherits(x, "foodnet_crm"))
    growth <- match.arg(growth)
    consumed <- x$consumed
    consumed[is.na(consumed)] <- 0
    total <- rowSums(consumed)
    squares <- rowSums(consumed^2)
    rates <- crm_growth(x, growth)
    dx <- x$biomass_change
    lacking <- names(rates)[is.na(rates) | is.na(dx) | dx <= 0 | total <= 0]
    if (length(lacking)) {
        stop_foodnet("a miaSim simulation needs a growth rate, a positive biomass change and a measured uptake ",
                     "for every taxon; ", paste(lacking, collapse = ", "), " lack one (see print(x)). Leave them ",
                     "out with crm_subset(x, taxa = ), or use crm_efficiency(x, scale = \"shares\") for a model ",
                     "that reads E as shares.")
    }
    out <- dx * total / (rates * squares)
    names(out) <- x$taxa
    out
}

#' A simulation's abundances in each taxon's growth curve unit
#'
#' @param x CRM parameters from [foodnet_listen()] or [foodnet_crm()].
#' @param abundance Abundances from a simulation built with [as_miasim()]: a vector with one value per
#'   taxon, or a matrix with taxa as rows (as `SummarizedExperiment::assay()` of miaSim's result).
#' @param growth The `growth` given to [as_miasim()].
#' @return The same, multiplied by [crm_scale()].
#' @export
crm_unscale <- function(x, abundance, growth = c("measured", "phase_floor")) {
    scale <- crm_scale(x, match.arg(growth))
    if (is.matrix(abundance)) {
        if (nrow(abundance) != length(scale)) {
            stop_foodnet("abundance needs one row per taxon (", length(scale), "), as assay() of miaSim's result; ",
                         "transpose it if taxa are its columns")
        }
        return(abundance * scale)
    }
    if (length(abundance) != length(scale)) stop_foodnet("abundance needs one value per taxon (", length(scale), ")")
    abundance * scale
}

# The growth rates a simulation uses, and a warning for any below its phase's own mean rate.
#' @noRd
crm_growth <- function(x, growth = "measured") {
    rates <- x$growth_rates
    phase <- x$phase_growth_rates
    if (is.null(phase)) return(rates)
    slow <- !is.na(phase) & !is.na(rates) & phase > 1.1 * rates
    if (growth == "phase_floor") return(ifelse(slow, phase, rates))
    if (any(slow)) {
        warning("growth rate below the phase's own mean rate for ", paste(names(rates)[slow], collapse = ", "),
                " (a fitted window took in the plateau); such a taxon grows too slowly in a simulation. ",
                "growth = \"phase_floor\" uses the phase's mean rate there.", call. = FALSE)
    }
    rates
}

#' Keep some taxa (and resources)
#'
#' The same parameters for fewer taxa or resources: for leaving out a taxon a simulation cannot use (no
#' growth rate, no biomass change, no measured uptake), or a resource. Every field indexed by taxon or
#' resource is subset together.
#'
#' @param x CRM parameters from [foodnet_listen()] or [foodnet_crm()].
#' @param taxa Taxa to keep, by name or position; all by default.
#' @param resources Resources to keep, by name or position; all by default.
#' @return CRM parameters of class `foodnet_crm`.
#' @export
crm_subset <- function(x, taxa = x$taxa, resources = x$resources) {
    stopifnot(inherits(x, "foodnet_crm"))
    keep_t <- if (is.character(taxa)) match(taxa, x$taxa) else taxa
    keep_r <- if (is.character(resources)) match(resources, x$resources) else resources
    if (anyNA(keep_t) || anyNA(keep_r)) stop_foodnet("unknown taxon or resource")
    for (name in c("consumed", "produced", "evidence_consumed", "evidence_produced")) {
        x[[name]] <- x[[name]][keep_t, keep_r, drop = FALSE]
    }
    for (name in c("growth_rates", "biomass_change", "biomass_start", "biomass_unit", "phase_hours",
                   "phase_growth_rates")) {
        x[[name]] <- x[[name]][keep_t]
    }
    x$initial <- x$initial[keep_r]
    kept <- x$taxa[keep_t]
    removed <- setdiff(x$taxa, kept)
    x$taxa <- kept
    x$resources <- x$resources[keep_r]
    p <- x$caveats$presence_only
    x$caveats$presence_only <- p[p$taxon %in% kept & p$resource %in% x$resources, , drop = FALSE]
    x$caveats$without_a_rate <- intersect(x$caveats$without_a_rate, kept)
    x$caveats$whole_run <- intersect(x$caveats$whole_run, kept)
    # the page's warnings were written for the whole search: keep those that name no taxon left out, and the
    # ones that name none (counts over the whole search are marked as such)
    if (length(x$caveats$warnings) && length(removed)) {
        names_any <- function(w) any(vapply(removed, function(t) grepl(t, w, fixed = TRUE), logical(1)))
        x$caveats$warnings <- c(paste("(for the whole search)", x$caveats$warnings[!vapply(x$caveats$warnings,
                                                                                           names_any, logical(1))]))
    }
    for (name in c("cautions", "inconclusive")) {
        rows <- x$caveats[[name]]
        if (NROW(rows)) x$caveats[[name]] <- rows[rows$taxon %in% kept & rows$resource %in% x$resources, , drop = FALSE]
    }
    x
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
#' starting abundances and the Monod constants. With the default `E`, the starting abundances are put into
#' each taxon's simulation unit ([crm_scale()]), and [crm_unscale()] turns the simulated abundances back. miaSim is not needed to call this, and no other simulator
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
#'   unit of that taxon's growth curve (`x$biomass_unit`); `x$biomass_start` is each taxon's own start. With
#'   your own `E`, they are passed as given.
#' @param monod_constant Monod constants in mM: one number for all, or a taxa by resources matrix.
#' @param E An efficiency matrix, by default [crm_efficiency()] of `x`.
#' @param missing_rate A growth rate for the taxa that have none, or `NULL` to stop when any does.
#' @param missing_resource A starting concentration for the resources that have none, or `NULL` to stop.
#' @param allow Caveats to accept: `"mixed_media"`, `"stationary_phase"`.
#' @param na What [crm_efficiency()] does with NA cells: `"stop"` (the default) or `"zero"`.
#' @param growth `"measured"` or `"phase_floor"`: the growth rates, as in [crm_scale()].
#' @return A list with `n_species`, `n_resources`, `names_species`, `names_resources`, `E`, `x0`,
#'   `resources`, `growth_rates`, `monod_constant` and `migration_p = 0`: miaSim adds random immigration
#'   even with `stochastic = FALSE` (miaSim 1.18's `perturb`), which in these units would swamp growth.
#' @examples
#' \dontrun{
#' args <- as_miasim(crm, x0 = crm$biomass_start, monod_constant = 1, missing_resource = 0, na = "zero")
#' set.seed(1)
#' tse <- do.call(miaSim::simulateConsumerResource, c(args, list(t_end = 48, t_store = 480)))
#' cells <- crm_unscale(crm, SummarizedExperiment::assay(tse))
#' }
#' @export
as_miasim <- function(x, x0, monod_constant, E = NULL, missing_rate = NULL, missing_resource = NULL,
                      allow = character(), na = c("stop", "zero"), growth = c("measured", "phase_floor")) {
    na <- match.arg(na)
    growth <- match.arg(growth)
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
    if (growth == "phase_floor") {
        # the same rule as the scale and E (crm_growth), and never a stand-in for a missing rate
        floored <- suppressWarnings(crm_growth(x, "phase_floor"))
        rates <- ifelse(is.na(x$growth_rates), rates, floored)
    }
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
    if (x$caveats$incomplete) {
        warning("these parameters are incomplete: records could not be read from mGrowthDB (print(x) lists ",
                "them)", call. = FALSE)
    }
    if (is.null(E)) {
        E <- crm_efficiency(x, na = na, growth = growth)
        x0 <- as.numeric(x0) / crm_scale(x, growth)       # into each taxon's own unit (crm_scale)
    }
    # migration_p = 0: miaSim adds random immigration even when stochastic is FALSE (its perturb() does not
    # scale that term by stochastic), and in these units one event is as large as a starting population
    list(n_species = n, n_resources = m, names_species = x$taxa, names_resources = x$resources, E = E,
         x0 = unname(as.numeric(x0)), resources = unname(resources), growth_rates = unname(rates),
         monod_constant = unname(K), migration_p = 0)
}

#' Simulate each taxon alone and compare with its monoculture
#'
#' Runs miaSim's consumer-resource model for each taxon on its own, from its measured starting abundance
#' and the medium's starting concentrations, over the length of its phase, and compares what it gains and
#' what it takes up and makes with what its monoculture did, and how long it takes to gain its measured
#' biomass with the length of its phase. With [crm_efficiency()]'s default scale the biomass and the
#' by-products follow the uptake by construction, so what this tests is the Monod constants: how much of
#' each resource is taken up, and how fast.
#'
#' @param x CRM parameters from [foodnet_listen()] or [foodnet_crm()].
#' @param monod_constant Monod constants in mM: one number, or a taxa by resources matrix.
#' @param missing_resource A starting concentration for resources with none (default 0).
#' @param na What [crm_efficiency()] does with NA cells: `"stop"` (the default) or `"zero"`.
#' @param growth `"measured"` or `"phase_floor"`: the growth rates, as in [crm_scale()].
#' @return A data frame: taxon, what (`"biomass"`, `"hours to grow"`, or a resource), direction, measured,
#'   simulated and their ratio. Biomass is in each growth curve's unit.
#' @examples
#' \dontrun{
#' crm_backcheck(crm, monod_constant = 1, na = "zero")
#' }
#' @export
crm_backcheck <- function(x, monod_constant, missing_resource = 0, na = c("stop", "zero"),
                          growth = c("measured", "phase_floor")) {
    na <- match.arg(na)
    growth <- match.arg(growth)
    stopifnot(inherits(x, "foodnet_crm"))
    if (!requireNamespace("miaSim", quietly = TRUE)) {
        stop_foodnet("crm_backcheck() runs miaSim: install it with BiocManager::install(\"miaSim\")")
    }
    if (missing(monod_constant)) stop_foodnet("give monod_constant, in mM; foodnet measures none")
    E <- crm_efficiency(x, na = na, growth = growth)
    scale <- crm_scale(x, growth)
    mu <- suppressWarnings(crm_growth(x, growth))
    n <- length(x$taxa)
    m <- length(x$resources)
    K <- if (length(monod_constant) == 1) matrix(monod_constant, n, m) else as.matrix(monod_constant)
    start <- crm_resources(x, missing = missing_resource)
    consumed <- x$consumed
    produced <- x$produced
    rows <- list()
    add <- function(...) rows[[length(rows) + 1]] <<- data.frame(..., stringsAsFactors = FALSE)
    for (i in seq_len(n)) {
        taxon <- x$taxa[i]
        hours <- x$phase_hours[[taxon]]
        x0 <- x$biomass_start[[taxon]]
        if (is.na(hours) || is.na(x0)) next
        # three phase lengths, so a taxon that grows too slowly shows how slowly
        steps <- max(10L, round(3 * hours / 0.1))
        tse <- miaSim::simulateConsumerResource(
            n_species = 1, n_resources = m, names_species = taxon, names_resources = x$resources,
            E = E[i, , drop = FALSE], x0 = x0 / scale[[i]], resources = unname(start),
            growth_rates = mu[[taxon]], monod_constant = K[i, , drop = FALSE],
            t_end = (steps + 1) * 0.1, t_store = steps + 1, migration_p = 0)
        abundance <- SummarizedExperiment::assay(tse)[1, ] * scale[[i]]
        times <- SummarizedExperiment::colData(tse)$time
        at <- which.min(abs(times - hours))
        grown <- abundance - abundance[1]
        reached <- which(grown >= x$biomass_change[[taxon]])
        add(taxon = taxon, what = "biomass", direction = "grown", measured = x$biomass_change[[taxon]],
            simulated = grown[at])
        add(taxon = taxon, what = "hours to grow", direction = "grown", measured = hours,
            simulated = if (length(reached)) times[reached[1]] else NA_real_)
        change <- S4Vectors::metadata(tse)$resources[at, seq_len(m)] - start
        for (j in seq_len(m)) {
            for (direction in c("consumed", "produced")) {
                measured <- if (direction == "consumed") consumed[i, j] else produced[i, j]
                simulated <- if (direction == "consumed") max(0, -change[j]) else max(0, change[j])
                if ((is.na(measured) || measured == 0) && simulated < 1e-6) next
                add(taxon = taxon, what = x$resources[j], direction = direction, measured = measured,
                    simulated = simulated)
            }
        }
    }
    out <- do.call(rbind, rows)
    if (is.null(out)) stop_foodnet("no taxon has a starting abundance and a phase length")
    out$ratio <- out$simulated / out$measured
    rownames(out) <- NULL
    out
}
