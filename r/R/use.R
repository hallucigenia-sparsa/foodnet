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
    need_crm(x)
    x$consumed
}

#' @rdname crm_consumed
#' @export
crm_produced <- function(x) {
    need_crm(x)
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
    need_crm(x)
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
    need_crm(x)
    initial <- x$initial
    if (!is.null(missing)) initial[is.na(initial)] <- missing
    initial
}

#' The efficiency matrix of a consumer-resource model
#'
#' A consumer-resource model such as miaSim's `simulateConsumerResource` takes an efficiency matrix `E`,
#' taxa by resources, positive where a taxon consumes a resource and negative where it produces one. This
#' builds one from the measured amounts, the consumed and produced matrices, which are never negative. foodnet's
#' single (signed) matrix has the opposite signs, production positive and consumption negative, as a change in
#' the medium reads; it is not an `E` and is never used here.
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
#' `scale = "none"` is consumed minus produced, in mM: foodnet's signed matrix with its signs turned to
#' miaSim's.
#'
#' `NA` cells have no size, and foodnet never writes one as zero, so by default (`na = "stop"`) this refuses
#' until you fill them yourself or say how. With `na = "zero"` they count as 0, and a warning says how many
#' of each kind were set to 0: never assayed, seen only in another medium (a real link), inconclusive
#' (measured, but its replicates or experiments disagree), without a phase, and, after [crm_phase()], a
#' second-window compound (`second_window`: its change is given once, in the exponential phase).
#'
#' With the default `scale = "miasim"`, like [as_miasim()], this refuses parameters that pool every medium or
#' come from the stationary phase, unless `allow` names them.
#'
#' @param x CRM parameters from [foodnet_listen()] or [foodnet_crm()].
#' @param scale `"miasim"`, `"shares"` or `"none"`, as above.
#' @param na `"stop"` or `"zero"`, as above.
#' @param allow Caveats to accept: `"mixed_media"`, `"stationary_phase"`.
#' @param growth `"measured"` (the default): the growth rates as foodnet fitted them. `"phase_floor"`: at
#'   least the mean rate each taxon's own curves show over the phase (see [crm_scale()]).
#' @return A numeric matrix, taxa by resources.
#' @examples
#' \dontrun{
#' E <- crm_efficiency(crm, na = "zero")
#' }
#' @export
crm_efficiency <- function(x, scale = c("miasim", "shares", "none"), na = c("stop", "zero"),
                           growth = c("measured", "phase_floor"), allow = character()) {
    need_crm(x)
    scale <- match.arg(scale)
    # miaSim's E reads every uptake as growth in one medium; "shares" and "none" are plain arithmetic on the
    # amounts for other models, as in 0.1.0
    if (scale == "miasim") refuse_caveats(x, allow)
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

# A clear stop when x is not the parameters themselves, e.g. crm_efficiency()'s matrix passed on (Karoline's
# test, 2026-10-07: crm_backcheck(crm_efficiency(crm, na = "zero"), ...)).
#' @noRd
need_crm <- function(x) {
    if (!inherits(x, "foodnet_crm")) {
        called <- tryCatch(deparse(sys.call(-1)[[1]]), error = function(e) "")
        # through lapply() or do.call() the caller has no name of ours: say it in general (a review)
        if (!called %in% getNamespaceExports("foodnet")) called <- "this function"
        stop_foodnet(called, "() takes the CRM parameters themselves (crm, from foodnet_listen() or ",
                     "foodnet_crm()), not ", if (is.matrix(x)) "a matrix such as crm_efficiency()'s E" else class(x)[1],
                     ": it builds what it needs from them, as in ", called, "(crm, ...)")
    }
    invisible(TRUE)
}

# A stop for parameters a CRM cannot be built from as they are, unless `allow` names the caveat.
#' @noRd
refuse_caveats <- function(x, allow = character()) {
    if (isTRUE(x$caveats$mixed_media) && !"mixed_media" %in% allow) {
        stop_foodnet("these values pool every medium (Ignore media differences): their starting concentrations ",
                     "describe no real medium. Run the search with one medium, or pass allow = \"mixed_media\".")
    }
    if (isTRUE(x$caveats$stationary_phase) && !"stationary_phase" %in% allow) {
        stop_foodnet("these are amounts from after the end of exponential growth, where cells may still grow, ",
                     "stop or die, while a CRM reads every uptake as growth. Use the exponential phase, or pass ",
                     "allow = \"stationary_phase\".")
    }
    invisible(TRUE)
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
                         "only in another medium, inconclusive, without a phase, or a second-window compound ",
                         "given in the other phase; crm$evidence_consumed and crm$evidence_produced say which). ",
                         "Fill them, or pass na = \"zero\" to the function you called to count them as 0.")
        }
        zeroed <- c(x$evidence_consumed[is.na(consumed)], x$evidence_produced[is.na(produced)])
        kinds <- table(zeroed)
        warning("NA cells counted as 0: ", paste0(kinds, " ", names(kinds), collapse = ", "),
                if ("presence_only" %in% names(kinds)) " (presence_only links are real; their size is unknown)",
                if ("second_window" %in% names(kinds)) " (second_window: measured over its own window, given once in the exponential phase)",
                call. = FALSE)
        consumed[is.na(consumed)] <- 0
        produced[is.na(produced)] <- 0
    }
    moving <- x$caveats$cautions
    moving <- if (is.null(moving)) character(0) else
        paste(moving$taxon, moving$resource)[grepl("still_changing", moving$cautions)]
    if (length(moving) && !isTRUE(x$caveats$stationary_phase)) {
        # a script that never prints the parameters still hears it; a stationary value holds that change, so it
        # is not underrated there (a review)
        warning(length(moving), " value(s) miss use that went on after the growth rate fell (still_changing), so E ",
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
    need_crm(x)
    growth <- match.arg(growth)
    consumed <- x$consumed
    consumed[is.na(consumed)] <- 0
    total <- rowSums(consumed)
    squares <- rowSums(consumed^2)
    rates <- crm_growth(x, growth)
    dx <- x$biomass_change
    lacking <- names(rates)[is.na(rates) | is.na(dx) | dx <= 0 | total <= 0]
    if (length(lacking) && isTRUE(x$caveats$stationary_phase)) {
        # after the end of exponential growth the biomass of most taxa falls, which no growth can scale; leaving
        # them out would leave a community of one (a review)
        stop_foodnet("miaSim's scale needs each taxon to gain biomass, and after the end of exponential growth ",
                     paste(lacking, collapse = ", "), if (length(lacking) == 1) " does" else " do", " not (or ",
                     "lack a growth rate or uptake). Stationary amounts cannot be scaled this way: use the ",
                     "exponential phase, or crm_efficiency(crm, scale = \"shares\", allow = \"stationary_phase\") ",
                     "for a model that reads E as shares.")
    }
    if (length(lacking)) {
        # name what each lacks, and give the line that leaves them out (Karoline's test, 2026-10-07: "see print(x)"
        # sent her looking)
        what <- vapply(lacking, function(t) {
            missing <- c(if (is.na(rates[[t]])) "growth rate",
                         if (is.na(dx[[t]])) "biomass change" else if (dx[[t]] <= 0) "biomass gain",
                         if (total[[t]] <= 0) "measured uptake")
            paste0(t, " (no ", paste(missing, collapse = ", no "), ")")
        }, character(1))
        stop_foodnet("a miaSim simulation needs a growth rate, a biomass gain and a measured uptake for every ",
                     "taxon, and ", paste(what, collapse = "; "), if (length(lacking) == 1) " lacks" else " lack",
                     " them. Leave ", if (length(lacking) == 1) "it" else "them", " out (crm being your parameters): ",
                     "crm <- crm_subset(crm, taxa = setdiff(crm$taxa, c(", paste0("\"", lacking, "\"", collapse = ", "),
                     "))); or use crm_efficiency(crm, scale = \"shares\") for a model that reads E as shares.")
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
#' @param growth The `growth` given to [as_miasim()]; with `args`, read from them.
#' @param args The list [as_miasim()] returned, so the unit is the one the simulation used (recommended):
#'   given `"phase_floor"` there and not here, a taxon whose rate it raised came back off by the ratio of
#'   the two rates.
#' @return The same, multiplied by [crm_scale()].
#' @export
crm_unscale <- function(x, abundance, growth = c("measured", "phase_floor"), args = NULL) {
    if (!is.null(args)) {
        used <- attr(args, "foodnet_growth")
        if (is.null(used)) stop_foodnet("args is not the list as_miasim() returned")
        if (is.na(used)) stop_foodnet("as_miasim() was given your own E, so its abundances are in your units already")
        if (!missing(growth) && !identical(match.arg(growth), used)) {
            stop_foodnet("growth = \"", match.arg(growth), "\", but the simulation used \"", used, "\"")
        }
        growth <- used
    }
    scale <- crm_scale(x, match.arg(growth, c("measured", "phase_floor")))
    made <- if (is.null(args)) NULL else attr(args, "foodnet_scale")
    if (!is.null(made) && (!identical(names(made), x$taxa) || !isTRUE(all.equal(unname(scale), unname(made))))) {
        stop_foodnet("args were made from other CRM parameters than these (their taxa or units differ): pass ",
                     "the crm you gave to as_miasim()")
    }
    if (is.matrix(abundance) && !is.null(rownames(abundance))) {
        # rows by name where they have names, so reordered rows keep their own unit (a review)
        if (!setequal(rownames(abundance), x$taxa)) stop_foodnet("the rows of abundance ", taxa_mismatch(rownames(abundance), x$taxa))
        return(abundance * scale[rownames(abundance)])
    }
    if (is.matrix(abundance)) {
        if (nrow(abundance) != length(scale)) {
            stop_foodnet("abundance needs one row per taxon (", length(scale), "), as assay() of miaSim's result; ",
                         "transpose it if taxa are its columns")
        }
        return(abundance * scale)
    }
    if (!is.null(names(abundance))) {
        # a named vector by name too, as one time point of the assay is (a review)
        if (!setequal(names(abundance), x$taxa)) stop_foodnet("the names of abundance ", taxa_mismatch(names(abundance), x$taxa))
        return(abundance * scale[names(abundance)])
    }
    if (length(abundance) != length(scale)) stop_foodnet("abundance needs one value per taxon (", length(scale), ")")
    abundance * scale
}

# What a set of names lacks or adds against the parameters' taxa, in words.
#' @noRd
taxa_mismatch <- function(given, taxa) {
    extra <- setdiff(given, taxa)
    lacking <- setdiff(taxa, given)
    paste0(if (length(extra)) paste0("name taxa these parameters lack (", paste(utils::head(extra, 3), collapse = ", "),
                                     ")") else "",
           if (length(extra) && length(lacking)) " and " else "",
           if (length(lacking)) paste0("leave out ", paste(utils::head(lacking, 3), collapse = ", ")) else "",
           "; give one value per taxon of crm$taxa")
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
    need_crm(x)
    keep_t <- if (is.character(taxa)) match(taxa, x$taxa) else taxa
    keep_r <- if (is.character(resources)) match(resources, x$resources) else resources
    if (anyNA(keep_t) || anyNA(keep_r)) stop_foodnet("unknown taxon or resource")
    # what a left-out resource was consumed or produced, so crm_electron_balance() withholds a share that would
    # miss it (a second-window compound is never counted, so leaving it out changes nothing)
    drop_r <- setdiff(seq_along(x$resources), keep_r)
    drop_r <- drop_r[!x$resources[drop_r] %in% x$second_window_resources]
    # (NA: no number for the taxon, though the medium held it when the phase began, so the balance still names it)
    limits <- if (length(x$detection_limits)) x$detection_limits else rep(x$detection_limit, length(x$resources))
    dropped <- function(b) {
        amount <- function(m) ifelse(is.na(m), 0, abs(m))
        gone <- amount(b$consumed[keep_t, drop_r, drop = FALSE]) + amount(b$produced[keep_t, drop_r, drop = FALSE])
        start <- if (is.null(b$phase_start)) matrix(NA_real_, length(x$taxa), length(x$resources)) else b$phase_start
        start <- ifelse(is.na(start), matrix(x$initial, length(x$taxa), length(x$resources), byrow = TRUE), start)
        held <- !is.na(start) & start > matrix(limits, length(x$taxa), length(x$resources), byrow = TRUE)
        unmeasured <- (is.na(b$consumed) & is.na(b$produced))[keep_t, drop_r, drop = FALSE]
        none <- unmeasured & held[keep_t, drop_r, drop = FALSE]
        gone[none] <- NA
        # no number and not held above the limit: it counted for nothing in the share (a review)
        gone[unmeasured & !none] <- -1
        if (!is.null(b$dropped)) gone <- cbind(b$dropped[keep_t, , drop = FALSE], gone)
        gone
    }
    # what such a resource could still hide in the floor: what the medium held below the limit, or NA where
    # nobody measured it in the medium (taken as absent, and named) (a review: dropping it raised the floor)
    gamma_all <- crm_chemistry(x)$degree_of_reduction
    below <- function(b) {
        start <- if (is.null(b$phase_start)) matrix(NA_real_, length(x$taxa), length(x$resources)) else b$phase_start
        start <- ifelse(is.na(start), matrix(x$initial, length(x$taxa), length(x$resources), byrow = TRUE), start)
        lim <- matrix(limits, length(x$taxa), length(x$resources), byrow = TRUE)
        g <- matrix(ifelse(!is.na(gamma_all) & gamma_all > 0, gamma_all, 0), length(x$taxa), length(x$resources),
                    byrow = TRUE)
        unmeasured <- is.na(b$consumed) & is.na(b$produced)
        held <- !is.na(start) & start > lim
        out <- ifelse(unmeasured & !held, pmin(start, lim) * g, 0)
        out <- out[keep_t, drop_r, drop = FALSE]
        if (!is.null(b$dropped_below)) out <- cbind(b$dropped_below[keep_t, , drop = FALSE], out)
        out
    }
    cultures_kept <- function(cm) if (is.null(cm)) NULL else lapply(cm[keep_t], function(m) m[, keep_r, drop = FALSE])
    x["dropped"] <- list(dropped(x))
    x["dropped_below"] <- list(below(x))
    # what each left-out resource could hide below its detection limit, in electrons (limit x degree of
    # reduction), for a balance that names it as incomplete
    gamma <- crm_chemistry(x)$degree_of_reduction
    hide <- ifelse(!is.na(gamma) & gamma > 0, limits * gamma, 0)[drop_r]
    names(hide) <- x$resources[drop_r]
    x$dropped_hidden <- c(x$dropped_hidden, hide)
    x["culture_changes"] <- list(cultures_kept(x$culture_changes))
    for (name in PHASE_MATRICES) {
        if (!is.null(x[[name]])) x[[name]] <- x[[name]][keep_t, keep_r, drop = FALSE]
    }
    for (name in c("growth_rates", PHASE_VECTORS)) {
        x[[name]] <- x[[name]][keep_t]
    }
    for (ph in names(x$other_phases)) {
        b <- x$other_phases[[ph]]
        b["dropped"] <- list(dropped(b))
        b["dropped_below"] <- list(below(b))
        b["culture_changes"] <- list(cultures_kept(b$culture_changes))
        for (name in PHASE_MATRICES) b[[name]] <- b[[name]][keep_t, keep_r, drop = FALSE]
        for (name in PHASE_VECTORS) b[[name]] <- b[[name]][keep_t]
        keep_rows <- function(rows) {
            if (NROW(rows)) rows[rows$taxon %in% x$taxa[keep_t] & rows$resource %in% x$resources[keep_r], , drop = FALSE]
            else rows
        }
        b$cautions <- keep_rows(b$cautions)
        b$inconclusive <- keep_rows(b$inconclusive)
        b$presence_only <- keep_rows(b$presence_only)
        b$whole_run <- intersect(b$whole_run, x$taxa[keep_t])
        b$biomass_falls <- intersect(b$biomass_falls, x$taxa[keep_t])
        x$other_phases[[ph]] <- b
    }
    x$initial <- x$initial[keep_r]
    if (length(x$detection_limits)) x$detection_limits <- x$detection_limits[keep_r]
    if (!is.null(x$chemistry)) x$chemistry <- x$chemistry[keep_r, , drop = FALSE]
    # one left out still qualifies the floor, which never counted it (a review)
    x$second_window_left_out <- union(x$second_window_left_out,
                                      setdiff(x$second_window_resources, x$resources[keep_r]))
    x$second_window_resources <- intersect(x$second_window_resources, x$resources[keep_r])
    if (!is.null(x$resource_phases)) x$resource_phases <- x$resource_phases[keep_r]
    kept <- x$taxa[keep_t]
    removed <- setdiff(x$taxa, kept)
    x$taxa <- kept
    x$resources <- x$resources[keep_r]
    p <- x$caveats$presence_only
    x$caveats$presence_only <- p[p$taxon %in% kept & p$resource %in% x$resources, , drop = FALSE]
    x$caveats$without_a_rate <- intersect(x$caveats$without_a_rate, kept)
    x$caveats$whole_run <- intersect(x$caveats$whole_run, kept)
    x$caveats$biomass_falls <- intersect(x$caveats$biomass_falls, kept)
    # the page's warnings were written for the whole search: keep those that name no taxon left out, and the
    # ones that name none (counts over the whole search are marked as such)
    if (length(x$caveats$warnings) && length(removed)) {
        # a removed taxon is named where its name appears outside every occurrence of a kept taxon's name
        # ("Escherichia coli" removed, "Escherichia coli LF82" kept; "Clostridium" removed, "[Clostridium]
        # scindens" kept)
        whole_name <- function(t, w) {
            at <- gregexpr(t, w, fixed = TRUE)[[1]]
            if (at[1] < 0) return(FALSE)
            covered <- integer(0)
            for (k in kept[nchar(kept) > nchar(t) & grepl(t, kept, fixed = TRUE)]) {
                ks <- gregexpr(k, w, fixed = TRUE)[[1]]
                if (ks[1] > 0) for (q in ks) covered <- c(covered, q:(q + nchar(k) - 1))
            }
            any(vapply(at, function(p) !all(p:(p + nchar(t) - 1) %in% covered), logical(1)))
        }
        names_any <- function(w) any(vapply(removed, whole_name, logical(1), w = w))
        left <- x$caveats$warnings[!vapply(x$caveats$warnings, names_any, logical(1))]
        mark <- "(for the whole search) "
        left <- ifelse(startsWith(left, mark), left, paste0(mark, left))
        x$caveats$warnings <- if (length(left)) left else character(0)
    }
    for (name in c("cautions", "inconclusive")) {
        rows <- x$caveats[[name]]
        if (NROW(rows)) x$caveats[[name]] <- rows[rows$taxon %in% kept & rows$resource %in% x$resources, , drop = FALSE]
    }
    x
}

# The fields that belong to one phase: what crm_phase() swaps and crm_subset() subsets, in the order a phase
# block holds them.
PHASE_MATRICES <- c("consumed", "produced", "evidence_consumed", "evidence_produced", "interval_start",
                    "interval_end", "consumed_lower", "consumed_upper", "produced_lower", "produced_upper",
                    "phase_start")
PHASE_VECTORS <- c("biomass_change", "biomass_start", "biomass_unit", "phase_hours", "phase_growth_rates")
PHASE_CAVEATS <- c("cautions", "inconclusive", "presence_only", "whole_run", "biomass_falls")

#' Switch to another phase
#'
#' A search for both phases sends the exponential phase, which a consumer-resource model describes, and the
#' stationary phase beside it: what changed after the end of exponential growth, to the last sample, which
#' the exponential phase misses (a compound taken up only late, for one). Cells may still grow there, stop,
#' or die (`x$biomass_change`, and the caveat `biomass_falls`). This returns the same parameters for that
#' phase: its amounts, evidence, intervals, bounds, growth and caveats; the taxa, resources, growth rates and
#' initial concentrations stay (the growth rates are the taxa's maximum rates, from exponential growth, in
#' every phase). A compound measured over the second time window is given once, in the exponential phase's
#' matrices, as its change over that window. [crm_efficiency()], [as_miasim()] and [crm_backcheck()] refuse stationary amounts
#' unless `allow = "stationary_phase"`, since a model reads every uptake as growth.
#'
#' @param x CRM parameters from [foodnet_listen()] or [foodnet_crm()].
#' @param phase The phase to switch to, one of `c(x$phase, names(x$other_phases))`.
#' @return CRM parameters of class `foodnet_crm`, for `phase`; the phase switched from is kept in
#'   `x$other_phases`, so switching back gives the original.
#' @export
crm_phase <- function(x, phase) {
    need_crm(x)
    if (identical(phase, x$phase)) return(x)
    b <- x$other_phases[[phase]]
    if (is.null(b)) {
        stop_foodnet("these parameters carry the ", paste(c(x$phase, names(x$other_phases)), collapse = " and "),
                     " phase only; search for both phases in foodnet to have the stationary one too")
    }
    current <- list(phase = x$phase)
    # the culture changes travel with their phase, and so does what crm_subset() left out, once it has
    for (name in intersect(c(PHASE_MATRICES, PHASE_VECTORS, "culture_changes", "dropped", "dropped_below"),
                           union(names(x), names(b)))) {
        current[name] <- list(x[[name]])
        x[name] <- list(b[[name]])
    }
    for (name in PHASE_CAVEATS) {
        current[[name]] <- x$caveats[[name]]
        x$caveats[name] <- list(b[[name]])
    }
    x$caveats$stationary_phase <- phase == "stationary"
    others <- x$other_phases
    others[[phase]] <- NULL
    others[[current$phase]] <- current
    x$other_phases <- others
    x$phase <- phase
    x
}

#' The README foodnet wrote with these numbers
#'
#' @param x CRM parameters from [foodnet_listen()] or [foodnet_crm()].
#' @param print Set to `FALSE` to return the text without printing it.
#' @return The text, invisibly when printed.
#' @export
crm_readme <- function(x, print = TRUE) {
    need_crm(x)
    if (print) {
        cat(x$readme)
        return(invisible(x$readme))
    }
    x$readme
}

#' Write the parameters as files
#'
#' The files of foodnet's download, for the phase `x` holds: `consumed.csv`, `produced.csv` and their
#' `evidence_*.csv`, `growth_rates.csv`, `initial_concentrations.csv`, `biomass.csv`, `intervals.csv` (the
#' hours each value was measured over) and `bounds.csv` (each amount's lowest and highest replicate), one row
#' per cell with a phase column, `chemistry.csv` and `electron_balance.csv` (from foodnet 0.3.0, [crm_chemistry()]
#' and [crm_electron_balance()], the latter with a phase column), and `README.txt`, which opens with the phase. After [crm_phase()] the
#' phase's files are named for it (`consumed_stationary.csv`, `intervals_stationary.csv`,
#' `README_stationary.txt`, ...), so both phases can be written to one folder; the growth rates and the
#' medium are the same in both. The tables hold what the parameters in R hold, so some have fewer columns than the download's
#' (no replicate counts, no rate sources, no value-medium header cell).
#'
#' @param x CRM parameters from [foodnet_listen()] or [foodnet_crm()].
#' @param dir Directory to write into. It is created when it does not exist.
#' @return The paths written, invisibly.
#' @export
crm_write <- function(x, dir) {
    need_crm(x)
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
    # the download names the CRM's own phase plainly and another phase by suffix: switched, x keeps the
    # exponential phase among the others
    own <- is.null(x$other_phases$exponential)
    suffix <- if (own) "" else paste0("_", x$phase)
    # after a switch, every file of the phase carries the suffix, so both phases can be written to one folder
    # (a review: the second write replaced the first's intervals, bounds and README); the growth rates and the
    # medium are the taxa's, whatever the phase
    phased <- c("consumed.csv", "produced.csv", "evidence_consumed.csv", "evidence_produced.csv", "biomass.csv",
                "intervals.csv", "bounds.csv", "electron_balance.csv")
    path <- function(name) file.path(dir, if (name %in% phased) sub("\\.csv$", paste0(suffix, ".csv"), name) else name)
    paths <- character(0)
    put <- function(table, name, row.names = FALSE) {
        p <- path(name)
        write.csv(table, p, row.names = row.names)
        paths <<- c(paths, p)
    }
    put(x$consumed, "consumed.csv", TRUE)
    put(x$produced, "produced.csv", TRUE)
    put(x$evidence_consumed, "evidence_consumed.csv", TRUE)
    put(x$evidence_produced, "evidence_produced.csv", TRUE)
    put(data.frame(taxon = x$taxa, growth_rate = unname(x$growth_rates), unit = x$growth_rate_unit,
                   stringsAsFactors = FALSE), "growth_rates.csv")
    put(data.frame(resource = x$resources, initial_mM = unname(x$initial), stringsAsFactors = FALSE),
        "initial_concentrations.csv")
    put(data.frame(taxon = x$taxa, biomass_start = unname(x$biomass_start),
                   biomass_change = unname(x$biomass_change), unit = unname(x$biomass_unit),
                   phase_hours = unname(x$phase_hours), stringsAsFactors = FALSE), "biomass.csv")
    cells <- expand.grid(i = seq_along(x$taxa), j = seq_along(x$resources))
    # a second-window compound was measured over that window, which the files say, as the download's do
    column_phase <- if (is.null(x$resource_phases)) rep(x$phase, length(x$resources)) else
        ifelse(x$resource_phases[x$resources] %in% "window", "window", x$phase)
    measured <- cells[!is.na(x$interval_start[as.matrix(cells)]), , drop = FALSE]
    put(data.frame(taxon = x$taxa[measured$i], resource = x$resources[measured$j],
                   phase = column_phase[measured$j],
                   start_h = x$interval_start[as.matrix(measured)], end_h = x$interval_end[as.matrix(measured)],
                   hours = x$interval_end[as.matrix(measured)] - x$interval_start[as.matrix(measured)],
                   stringsAsFactors = FALSE), "intervals.csv")
    bounds <- do.call(rbind, lapply(c("consumed", "produced"), function(d) {
        low <- x[[paste0(d, "_lower")]]
        kept <- cells[!is.na(low[as.matrix(cells)]), , drop = FALSE]
        at <- as.matrix(kept)
        data.frame(taxon = x$taxa[kept$i], resource = x$resources[kept$j], phase = column_phase[kept$j],
                   direction = rep(d, nrow(kept)),
                   value_mM = x[[d]][at], lower_mM = low[at], upper_mM = x[[paste0(d, "_upper")]][at],
                   stringsAsFactors = FALSE)
    }))
    put(bounds, "bounds.csv")
    if (!is.null(x$chemistry_source)) {
        # as foodnet's download: every resource, and the taxa with something counted or a share withheld
        put(crm_chemistry(x), "chemistry.csv")
        b <- crm_electron_balance(x)
        b <- b[!is.na(b$withheld) | (!is.na(b$consumed_electrons_mM) &
                                     (b$consumed_electrons_mM != 0 | b$produced_electrons_mM != 0)), , drop = FALSE]
        put(cbind(b[, 1, drop = FALSE], phase = rep(x$phase, nrow(b)), b[, -1, drop = FALSE]), "electron_balance.csv")
    }
    readme <- file.path(dir, if (own) "README.txt" else paste0("README", suffix, ".txt"))
    writeLines(c(sprintf("These files hold the %s phase%s.", x$phase,
                         if (own) "" else " (switched to with crm_phase(); the README below describes the search)"),
                 "", x$readme), readme)
    invisible(c(paths, readme))
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
#' real medium) or come from after the end of exponential growth (the stationary phase), unless `allow`
#' names them.
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
#'   even with `stochastic = FALSE` (miaSim 1.18's `perturb`), which in these units would swamp growth. Its
#'   other noise (drift, epochs, external events) is off unless `stochastic = TRUE`, and measurement noise
#'   unless `error_variance > 0`. To explore noise, change them in this list before the call, since it
#'   already holds `migration_p` (`args$migration_p <- 0.01; args$stochastic <- TRUE`); to turn it all off
#'   again, `args$migration_p <- 0; args$stochastic <- FALSE; args$error_variance <- 0`. Leave `norm = FALSE`: relative
#'   abundances cannot be turned back by [crm_unscale()].
#' @examples
#' \dontrun{
#' args <- as_miasim(crm, x0 = crm$biomass_start, monod_constant = 1, missing_resource = 0, na = "zero")
#' tse <- do.call(miaSim::simulateConsumerResource, c(args, list(t_end = 48, t_store = 480)))
#' cells <- crm_unscale(crm, SummarizedExperiment::assay(tse), args = args)
#' }
#' @export
as_miasim <- function(x, x0, monod_constant, E = NULL, missing_rate = NULL, missing_resource = NULL,
                      allow = character(), na = c("stop", "zero"), growth = c("measured", "phase_floor")) {
    na <- match.arg(na)
    growth <- match.arg(growth)
    need_crm(x)
    refuse_caveats(x, allow)
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
    # a taxon a simulation cannot use is named first, with what it lacks, rather than as an NA in x0 (a review:
    # "x0 needs one number per taxon" sent a user looking at x0)
    if (is.null(E)) suppressWarnings(crm_scale(x, growth))       # crm_efficiency() below says what it warns
    x0 <- if (!is.null(names(x0)) && all(x$taxa %in% names(x0))) x0[x$taxa] else x0
    if (length(x0) != n) stop_foodnet("x0 needs one number per taxon (", n, "), in the order of crm$taxa or named")
    if (anyNA(x0)) stop_foodnet("x0 is NA for ", paste(x$taxa[is.na(x0)], collapse = ", "), ": give a starting ",
                                "abundance, or leave the taxon out with crm_subset()")
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
    own_E <- !is.null(E)
    if (is.null(E)) {
        E <- crm_efficiency(x, na = na, growth = growth, allow = allow)
        x0 <- as.numeric(x0) / crm_scale(x, growth)       # into each taxon's own unit (crm_scale)
    }
    # migration_p = 0: miaSim adds random immigration even when stochastic is FALSE (its perturb() does not
    # scale that term by stochastic), and in these units one event is as large as a starting population
    args <- list(n_species = n, n_resources = m, names_species = x$taxa, names_resources = x$resources, E = E,
                 x0 = unname(as.numeric(x0)), resources = unname(resources), growth_rates = unname(rates),
                 monod_constant = unname(K), migration_p = 0)
    # the unit the abundances are in, for crm_unscale(args = ): with "phase_floor" given here and not there, a
    # taxon came back off by the ratio of its two rates (a review); NA with your own E (no unit of ours).
    # c(args, list(...)) drops it, so do.call() never passes it to miaSim.
    attr(args, "foodnet_growth") <- if (own_E) NA_character_ else growth
    # and the unit itself, so the args of other parameters cannot unscale this run unnoticed (a review)
    if (!own_E) attr(args, "foodnet_scale") <- suppressWarnings(crm_scale(x, growth))
    args
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
#' @param allow Caveats to accept, as in [crm_efficiency()].
#' @return A data frame: taxon, what (`"biomass"`, `"hours to grow"`, or a resource), direction, measured,
#'   simulated and their ratio. Biomass is in each growth curve's unit. "biomass" compares what the taxon gains
#'   over its phase; "hours to grow" compares the time it takes to gain its measured biomass with the length
#'   of its phase (NA: it never gains all of it within three phase lengths, either too slowly or levelling off
#'   just short; the biomass row says which). Both rows judge the Monod constants only over the exponential
#'   phase, whose length is the time each taxon grew. Over a time window neither does: the measured hours are
#'   the window's length, and while the constants lie well below the resource concentrations the biomass
#'   ratio is near 1 whatever they are, so this warns.
#' @examples
#' \dontrun{
#' crm_backcheck(crm, monod_constant = 1, na = "zero")
#' }
#' @export
crm_backcheck <- function(x, monod_constant, missing_resource = 0, na = c("stop", "zero"),
                          growth = c("measured", "phase_floor"), allow = character()) {
    na <- match.arg(na)
    growth <- match.arg(growth)
    need_crm(x)
    if (!requireNamespace("miaSim", quietly = TRUE)) {
        stop_foodnet("crm_backcheck() runs miaSim: install it with BiocManager::install(\"miaSim\")")
    }
    if (missing(monod_constant)) stop_foodnet("give monod_constant, in mM; foodnet measures none")
    if (identical(x$phase, "window")) {
        warning("these parameters are over a time window, whose length is not the time each taxon grew: the ",
                "hours-to-grow ratios do not judge the Monod constants, and neither do the biomass ratios while ",
                "the constants lie well below the resource concentrations. Search the same taxa with the ",
                "exponential phase to choose them.", call. = FALSE)
    }
    E <- crm_efficiency(x, na = na, growth = growth, allow = allow)
    scale <- suppressWarnings(crm_scale(x, growth))             # crm_efficiency() above has warned once
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

#' The resources' chemistry
#'
#' Each resource's formula and charge from ChEBI, by the ChEBI id foodnet keys it under, its carbon atoms and
#' its degree of reduction: the electrons per molecule relative to CO2, H2O, NH3, H2SO4, H3PO4 and H+,
#' `4C + H - 2O - 3N + 6S + 5P - charge` (Roels 1983). The charge term makes an acid and its conjugate base,
#' which foodnet joins, the same. A ChEBI class without a formula (succinate, fructose) takes the formula of a
#' named form (`formula_from`); a formula that cannot be evaluated is NA, and `note` says why.
#'
#' @param x CRM parameters from [foodnet_listen()] or [foodnet_crm()].
#' @return A data frame with a row per resource: `resource`, `chebi_id`, `formula`, `charge`, `carbon`,
#'   `degree_of_reduction`, `per_cmol` (per carbon atom), `formula_from` and `note`. All NA for parameters
#'   from foodnet before 0.3.0.
#' @export
crm_chemistry <- function(x) {
    need_crm(x)
    if (is.null(x$chemistry)) chemistry_frame(NULL, x$resources) else x$chemistry
}

#' Each taxon's electron balance
#'
#' The electrons in a taxon's measured by-products over the electrons in what it consumed of the measured
#' resources, `sum(produced * degree of reduction) / sum(consumed * degree of reduction)`, over the phase of
#' `x` ([crm_phase()] switches). It is a check and changes no value, and it sees one part of each side of the
#' balance: electrons also come from substrates nobody measured (peptides and amino acids of a rich medium,
#' dihydrogen or formate taken up by hydrogen-using taxa) and go to biomass and to by-products nobody
#' measured. So it has no expected side of 1: above 1, unmeasured substrates (or a measurement problem) gave
#' more than biomass and unmeasured products took; near 1 is no proof the balance closes; and 1 minus the
#' share is not a biomass yield.
#'
#' The share is withheld when a resource the taxon consumed or produced has no degree of reduction, or was left
#' out with [crm_subset()], since leaving out a substrate makes the share too high and a product too low.
#' It is also withheld when the electrons taken up are within the detection limits of the resources counted and
#' of those of the medium it has no number for (each limit times its degree of reduction, summed), since the
#' share would then be undetermined; `share_at_least` then holds what it is at least (electrons out over electrons
#' in plus the most the resources could hide: for the pooled amounts, a counted one its detection limit and one
#' without a number what the medium held when the phase began; for a culture that measured all of them, its own
#' changes, each hiding its limit; the lowest of these). And it is withheld when it lies outside its own cultures' range, which the
#' cultures left out of the range then drive. Compounds of
#' the second time window are not counted (measured over another window). The range is the lowest and highest
#' culture's own share over the same resources as the share, from the cultures that measured all of them and
#' took up more than those limits; a culture's change inside the detection limit is 0, as in the matrices. What
#' [crm_subset()] left out is kept in `x$dropped`.
#'
#' @param x CRM parameters from [foodnet_listen()] or [foodnet_crm()].
#' @return A data frame with a row per taxon: `consumed_electrons_mM` and `produced_electrons_mM` (mM of
#'   electrons; 0 where nothing was counted), `share`, `share_lower` and `share_upper` (the lowest and highest
#'   culture's share; NA with fewer than two such cultures, without them in parameters from before 0.3.0, or
#'   with an electron acceptor, a negative degree of reduction, counted), `cultures` (how many),
#'   `cultures_left_out` (cultures that measured those resources but took up no more electrons than the
#'   detection limits hide), `withheld` (why
#'   there is no share although there are numbers), `incomplete` (resources of the medium without a number for
#'   the taxon when the phase began, which the share leaves out) and `not_counted` (every resource left out, and
#'   why). All NA
#'   without degrees of reduction, or when the values are booleans.
#' @export
crm_electron_balance <- function(x) {
    need_crm(x)
    columns <- data.frame(taxon = character(0), consumed_electrons_mM = numeric(0),
                          produced_electrons_mM = numeric(0), share = numeric(0), share_lower = numeric(0),
                          share_upper = numeric(0), share_at_least = numeric(0), cultures = integer(0),
                          cultures_left_out = integer(0),
                          withheld = character(0),
                          incomplete = character(0), not_counted = character(0), stringsAsFactors = FALSE)
    if (!length(x$taxa)) return(columns)
    gamma <- crm_chemistry(x)$degree_of_reduction
    if (all(is.na(gamma)) || identical(x$values, "booleans")) {
        out <- columns[rep(NA_integer_, length(x$taxa)), , drop = FALSE]
        out$taxon <- x$taxa
        out$cultures <- 0L
        out$cultures_left_out <- 0L
        rownames(out) <- NULL
        return(out)
    }
    counted <- !x$resources %in% x$second_window_resources
    # a resource is of the medium when the taxon's cultures held it above the limit when the phase began, or,
    # where they did not measure it, the medium did at its first sample (Karoline, 2026-10-08)
    start <- if (is.null(x$phase_start)) matrix(NA_real_, length(x$taxa), length(x$resources)) else x$phase_start
    medium <- matrix(x$initial, nrow = length(x$taxa), ncol = length(x$resources), byrow = TRUE)
    start <- ifelse(is.na(start), medium, start)
    limits <- if (length(x$detection_limits)) x$detection_limits else rep(x$detection_limit, length(x$resources))
    in_medium_of <- function(i) !is.na(start[i, ]) & start[i, ] > limits
    rows <- lapply(seq_along(x$taxa), function(i) {
        e_in <- 0
        e_out <- 0
        used <- integer(0)
        skipped <- character(0)
        incomplete <- character(0)
        missing <- FALSE
        for (j in seq_along(x$resources)) {
            c_v <- x$consumed[i, j]
            p_v <- x$produced[i, j]
            if (!counted[j]) {
                skipped <- c(skipped, paste0(x$resources[j], ": measured over the second time window"))
                next
            }
            if (is.na(c_v) && is.na(p_v)) {
                skipped <- c(skipped, paste0(x$resources[j], ": no number"))
                if (in_medium_of(i)[j]) incomplete <- c(incomplete, x$resources[j])
                next
            }
            if (is.na(gamma[j])) {
                if (any(!is.na(c(c_v, p_v)) & c(c_v, p_v) != 0)) {
                    skipped <- c(skipped, paste0(x$resources[j], ": no degree of reduction"))
                    missing <- TRUE
                }
                next
            }
            if (is.na(c_v)) {
                skipped <- c(skipped, paste0(x$resources[j], ": consumed: no number"))
                if (in_medium_of(i)[j]) incomplete <- c(incomplete, x$resources[j])
            }
            if (is.na(p_v)) skipped <- c(skipped, paste0(x$resources[j], ": produced: no number"))
            e_in <- e_in + (if (is.na(c_v)) 0 else c_v) * gamma[j]
            e_out <- e_out + (if (is.na(p_v)) 0 else p_v) * gamma[j]
            used <- c(used, j)
        }
        withheld <- NA_character_
        if (missing) {
            withheld <- paste("no degree of reduction for a resource it consumed or produced, so the share would",
                              "leave it out (see not_counted)")
        }
        gone <- x$dropped
        missing_left <- if (is.null(gone)) character(0) else colnames(gone)[is.na(gone[i, ])]
        if (length(missing_left)) {
            # a resource left out that the medium held when the phase began, without a number for the taxon
            incomplete <- c(incomplete, paste0(missing_left, " (left out with crm_subset())"))
        }
        if (is.na(withheld) && !is.null(gone) && any(!is.na(gone[i, ]) & gone[i, ] > 0)) {
            withheld <- paste0("left out with crm_subset() although consumed or produced: ",
                               paste(colnames(gone)[!is.na(gone[i, ]) & gone[i, ] > 0], collapse = ", "))
        }
        # the uptake every counted resource could hide below its detection limit, in electrons (Karoline,
        # 2026-10-08: "Withhold, say why")
        # (every resource counted, and every one the medium held without a number for the taxon)
        hidden <- union(used, match(intersect(incomplete, x$resources), x$resources))
        # and those crm_subset() left out without a number or measured as 0, which still hide as much (reviews:
        # leaving them out brought back a share withheld for it)
        zero_left <- if (is.null(gone)) character(0) else colnames(gone)[!is.na(gone[i, ]) & gone[i, ] == 0]
        noise <- sum((limits * gamma)[hidden][!is.na(gamma[hidden]) & gamma[hidden] > 0]) +
            sum(x$dropped_hidden[c(missing_left, zero_left)], na.rm = TRUE)
        small <- FALSE
        floor <- NA_real_
        if (is.na(withheld) && e_in > 0 && e_in <= noise) {
            small <- TRUE
            withheld <- paste0("the electrons taken up (", sprintf("%.3g", e_in), " mM) are ",
                               "within the detection limits of the resources counted or missing (",
                               sprintf("%.3g", noise), " mM), so the share would be undetermined")
        }
        share <- if (e_in > 0 && is.na(withheld)) e_out / e_in else NA_real_
        # each culture's own share, over the resources the share counts, from those that measured all of them
        measured <- matrix(numeric(0), ncol = 2)
        cm <- x$culture_changes[[i]]
        if ((!is.na(share) || small) && !is.null(cm) && nrow(cm) && all(gamma[used] >= 0)) {
            for (k in seq_len(nrow(cm))) {
                v <- cm[k, used]
                if (anyNA(v)) next
                measured <- rbind(measured, c(sum(pmax(0, -v) * gamma[used]), sum(pmax(0, v) * gamma[used])))
            }
        }
        judged <- measured[, 1] > noise & measured[, 1] > 0
        shares <- measured[judged, 2] / measured[judged, 1]
        left_out <- as.integer(sum(!judged))
        ranged <- length(shares) >= 2
        low <- if (ranged) min(shares) else NA_real_
        high <- if (ranged) max(shares) else NA_real_
        if (!is.na(share) && ranged && !(share >= low - 1e-9 && share <= high + 1e-9)) {
            # outside its cultures' own range, the others drive it (Karoline, 2026-10-08)
            withheld <- paste0("the share (", sprintf("%.3g", share), ") lies outside its cultures' own (",
                               sprintf("%.3g", low), " to ", sprintf("%.3g", high),
                               "), so cultures outside the range (too little uptake to judge, or not every ",
                               "resource counted) drive it")
            share <- NA_real_
        }
        if (small) {
            # what the measured numbers still say (Karoline, 2026-10-08: "Report 'at least X'"): a counted
            # resource hides at most its detection limit, one without a number at most what the medium held when
            # the phase began; the lowest of the pooled amounts' and each culture's own
            inc <- match(intersect(incomplete, x$resources), x$resources)
            miss <- inc[!is.na(gamma[inc]) & gamma[inc] > 0]
            if (!length(missing_left) && !anyNA(gamma[inc]) && all(gamma[inc] >= 0) && !anyNA(start[i, miss])) {
                limit_e <- function(j) sum((limits * gamma)[j][gamma[j] > 0])
                # a resource without a number that the medium held below the limit hides at most what it held; one
                # nobody measured in the medium is taken as absent, and the reason says so (a review)
                unseen <- which(counted & is.na(x$consumed[i, ]) & is.na(x$produced[i, ]) &
                                    !x$resources %in% incomplete & !is.na(gamma) & gamma > 0)
                absent <- unseen[is.na(start[i, unseen])]
                seen <- setdiff(unseen, absent)
                # (named, also with one column: a single element would lose its name)
                gone_below <- if (is.null(x$dropped_below)) numeric(0) else
                    stats::setNames(x$dropped_below[i, ], colnames(x$dropped_below))
                below <- sum(pmin(start[i, seen], limits[seen]) * gamma[seen]) + sum(gone_below, na.rm = TRUE)
                most <- limit_e(used) + sum(start[i, miss] * gamma[miss]) + below +
                    sum(x$dropped_hidden[zero_left], na.rm = TRUE)
                floors <- e_out / (e_in + most)
                # a culture by its own numbers, where it measured every resource counted and every one the pooled
                # values lack, each hiding at most its limit (a review)
                every <- sort(union(used, miss))
                if (!is.null(cm) && nrow(cm)) {
                    for (k in seq_len(nrow(cm))) {
                        v <- cm[k, every]
                        if (anyNA(v)) next
                        floors <- c(floors, sum(pmax(0, v) * gamma[every]) /
                                        (sum(pmax(0, -v) * gamma[every]) + limit_e(every) + below +
                                             sum(x$dropped_hidden[zero_left], na.rm = TRUE)))
                    }
                }
                floor <- min(floors)
                notes <- c(if (any(!counted) || length(x$second_window_left_out))
                               "leaving out the second time window's compounds",
                           if (length(absent) || any(is.na(gone_below)))
                               paste0("taking ", paste(c(x$resources[absent], names(gone_below)[is.na(gone_below)]),
                                                       collapse = ", "),
                                                      ", which nobody measured in this medium, as absent"))
                withheld <- paste0(withheld, "; it is at least ", sprintf("%.3g", floor),
                                   if (length(notes)) paste0(" (", paste(notes, collapse = ", and "), ")") else "")
            } else {
                withheld <- paste0(withheld, "; a resource without a number has no bound, so no floor either")
            }
        }
        data.frame(taxon = x$taxa[i], consumed_electrons_mM = e_in, produced_electrons_mM = e_out, share = share,
                   share_lower = low, share_upper = high, share_at_least = floor,
                   cultures = length(shares), cultures_left_out = left_out, withheld = withheld,
                   incomplete = paste(incomplete, collapse = "; "),
                   not_counted = paste(skipped, collapse = "; "), stringsAsFactors = FALSE)
    })
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    out
}
