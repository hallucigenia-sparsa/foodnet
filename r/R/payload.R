# The CRM payload foodnet sends, and the object it becomes in R.
#
# As in grownet's package, the caveats travel as data beside the numbers, and this file is where they
# become part of the object rather than prose a reader may never open: `print` shows them every time, and
# `crm_efficiency` and `as_miasim` act on them.

CRM_FORMAT <- "foodnet.crm/v1"
# what this package still reads: v0 lacks the biomass changes, so crm_efficiency() then builds only the
# "shares" and "none" scales
CRM_FORMATS <- c("foodnet.crm/v0", "foodnet.crm/v1")

#' @noRd
stop_foodnet <- function(...) stop(paste0(...), call. = FALSE)

#' @noRd
as_number <- function(x) {
    if (is.null(x)) NA_real_ else as.numeric(x)
}

#' @noRd
chr <- function(x, default = "") {
    if (is.null(x) || !length(x)) default else as.character(x)[1]
}

#' @noRd
chr_vector <- function(x) {
    if (is.null(x) || !length(x)) character(0) else vapply(x, as.character, character(1))
}

# A list of rows (from JSON, never simplified, so a null stays a null) to a matrix with names.
#' @noRd
as_matrix <- function(rows, rownames, colnames, kind = "numeric") {
    n <- length(rownames)
    m <- length(colnames)
    convert <- if (kind == "numeric") as_number else function(v) chr(v, NA_character_)
    template <- if (kind == "numeric") numeric(1) else character(1)
    cells <- if (n && m) vapply(unlist(rows, recursive = FALSE), convert, template) else template[0]
    matrix(cells, nrow = n, ncol = m, byrow = TRUE, dimnames = list(rownames, colnames))
}

# A vector with names.
#' @noRd
named <- function(values, labels) {
    names(values) <- labels
    values
}

# A taxa by resources matrix of numbers that a payload from foodnet before 0.2.0 lacks (intervals in hours,
# bounds in mM): all NA then.
#' @noRd
hours_matrix <- function(rows, taxa, resources) {
    if (is.null(rows)) {
        return(matrix(NA_real_, nrow = length(taxa), ncol = length(resources), dimnames = list(taxa, resources)))
    }
    as_matrix(rows, taxa, resources)
}

# One number per taxon from a payload field, named by taxon (NA where the field is missing).
#' @noRd
taxon_numbers <- function(obj, field, taxa) {
    values <- if (length(obj[[field]])) vapply(obj[[field]], as_number, numeric(1)) else rep(NA_real_, length(taxa))
    names(values) <- taxa
    values
}

# One phase of the payload beside the CRM's own (foodnet's other_phases): the same fields, by the same
# taxa and resources.
#' @noRd
as_phase_block <- function(b, taxa, resources) {
    units <- if (length(b$biomass_unit)) vapply(b$biomass_unit, chr, character(1), NA_character_) else
        rep(NA_character_, length(taxa))
    names(units) <- taxa
    list(phase = chr(b$phase),
         consumed = as_matrix(b$consumed, taxa, resources),
         produced = as_matrix(b$produced, taxa, resources),
         evidence_consumed = as_matrix(b$evidence_consumed, taxa, resources, "character"),
         evidence_produced = as_matrix(b$evidence_produced, taxa, resources, "character"),
         interval_start = hours_matrix(b$interval_start_h, taxa, resources),
         interval_end = hours_matrix(b$interval_end_h, taxa, resources),
         consumed_lower = hours_matrix(b$consumed_lower, taxa, resources),
         consumed_upper = hours_matrix(b$consumed_upper, taxa, resources),
         produced_lower = hours_matrix(b$produced_lower, taxa, resources),
         produced_upper = hours_matrix(b$produced_upper, taxa, resources),
         biomass_change = taxon_numbers(b, "biomass_change", taxa),
         biomass_start = taxon_numbers(b, "biomass_start", taxa),
         biomass_unit = units,
         phase_hours = taxon_numbers(b, "phase_hours", taxa),
         phase_growth_rates = taxon_numbers(b, "phase_growth_rates", taxa),
         cautions = cautions_rows(b$cautions),
         inconclusive = caveats_rows(b$inconclusive),
         presence_only = presence_rows(b$presence_only),
         whole_run = chr_vector(b$whole_run),
         biomass_falls = chr_vector(b$biomass_falls))
}

# The links seen only in another medium, as a data frame (taxon, resource, direction, media).
#' @noRd
presence_rows <- function(presence) {
    if (!length(presence)) {
        return(data.frame(taxon = character(0), resource = character(0), direction = character(0),
                          media = character(0), stringsAsFactors = FALSE))
    }
    data.frame(taxon = vapply(presence, function(p) chr(p$taxon), character(1)),
               resource = vapply(presence, function(p) chr(p$resource), character(1)),
               direction = vapply(presence, function(p) chr(p$direction), character(1)),
               media = vapply(presence, function(p) paste(chr_vector(p$media), collapse = "; "), character(1)),
               stringsAsFactors = FALSE)
}

# A payload to the object this package works with.
#' @noRd
as_foodnet_crm <- function(payload) {
    if (!is.list(payload) || is.null(payload$format)) {
        stop_foodnet("this is not a foodnet CRM payload: it has no format field")
    }
    if (!chr(payload$format) %in% CRM_FORMATS) {
        # a newer format can change what a number means, so it is not read as if it were this one (a review: the
        # 0.1.0 package read 0.2.0's parameters with a warning and filled their inconclusive cells with 0)
        major <- suppressWarnings(as.integer(sub("^foodnet\\.crm/v([0-9]+).*$", "\\1", chr(payload$format))))
        newest <- max(as.integer(sub("^foodnet\\.crm/v", "", CRM_FORMATS)))
        if (!is.na(major) && major > newest) {
            stop_foodnet("these parameters are ", chr(payload$format), ", from a newer foodnet than this package ",
                         "reads (", paste(CRM_FORMATS, collapse = ", "), "): install the R package of that ",
                         "foodnet's version; its page prints the line")
        }
        warning("this payload says it is ", chr(payload$format), ", and this package reads ",
                paste(CRM_FORMATS, collapse = " and "), "; reading it anyway", call. = FALSE)
    }
    taxa <- chr_vector(payload$taxa)
    resources <- chr_vector(payload$resources)
    rates <- vapply(payload$growth_rates, as_number, numeric(1))
    names(rates) <- taxa
    initial <- vapply(payload$initial_concentrations, as_number, numeric(1))
    names(initial) <- resources
    per_taxon <- function(field) {
        values <- if (length(payload[[field]])) vapply(payload[[field]], as_number, numeric(1)) else
            rep(NA_real_, length(taxa))
        names(values) <- taxa
        values
    }
    units <- if (length(payload$biomass_unit)) vapply(payload$biomass_unit, chr, character(1), NA_character_) else
        rep(NA_character_, length(taxa))
    names(units) <- taxa
    inconclusive <- caveats_rows(payload$caveats$inconclusive)
    caveats <- payload$caveats
    presence <- presence_rows(caveats$presence_only)
    structure(
        list(taxa = taxa,
             resources = resources,
             consumed = as_matrix(payload$consumed, taxa, resources),
             produced = as_matrix(payload$produced, taxa, resources),
             evidence_consumed = as_matrix(payload$evidence_consumed, taxa, resources, "character"),
             evidence_produced = as_matrix(payload$evidence_produced, taxa, resources, "character"),
             interval_start = hours_matrix(payload$interval_start_h, taxa, resources),
             interval_end = hours_matrix(payload$interval_end_h, taxa, resources),
             consumed_lower = hours_matrix(payload$consumed_lower, taxa, resources),
             consumed_upper = hours_matrix(payload$consumed_upper, taxa, resources),
             produced_lower = hours_matrix(payload$produced_lower, taxa, resources),
             produced_upper = hours_matrix(payload$produced_upper, taxa, resources),
             bounds = chr(payload$bounds),
             resource_phases = named(if (length(payload$resource_phases))
                 vapply(payload$resource_phases, chr, character(1)) else rep(chr(payload$phase), length(resources)),
                 resources),
             growth_rates = rates,
             growth_rate_unit = chr(payload$growth_rate_unit, "1/h"),
             growth_rate_detail = payload$growth_rate_detail,
             initial = initial,
             biomass_change = per_taxon("biomass_change"),
             biomass_start = per_taxon("biomass_start"),
             biomass_unit = units,
             phase_hours = per_taxon("phase_hours"),
             phase_growth_rates = per_taxon("phase_growth_rates"),
             phase = chr(payload$phase),
             other_phases = lapply(payload$other_phases, as_phase_block, taxa = taxa, resources = resources),
             values = chr(payload$values, "mM"),
             detection_limit = as_number(payload$detection_limit_mM),
             caveats = list(presence_only = presence,
                            conflicts = chr_vector(caveats$conflicts),
                            duplicates = chr_vector(caveats$duplicates),
                            without_a_rate = chr_vector(caveats$without_a_rate),
                            media = chr_vector(caveats$media),
                            value_rule = chr(caveats$value_rule),
                            searched_both_phases = isTRUE(caveats$searched_both_phases),
                            mixed_media = isTRUE(caveats$mixed_media),
                            stationary_phase = isTRUE(caveats$stationary_phase),
                            whole_run = chr_vector(caveats$whole_run),
                            biomass_falls = chr_vector(caveats$biomass_falls),
                            caution_tiers = unlist(caveats$caution_tiers),
                            inconclusive = inconclusive,
                            incomplete = isTRUE(caveats$incomplete),
                            errors = chr_vector(caveats$errors),
                            warnings = chr_vector(caveats$warnings),
                            cautions = cautions_rows(caveats$cautions)),
             readme = chr(payload$readme),
             tool = chr(payload$tool, "foodnet"),
             tool_version = chr(payload$tool_version),
             derived_at = chr(payload$derived_at),
             source_db = chr(payload$source_db),
             studies = chr_vector(payload$studies),
             settings = payload$settings),
        class = "foodnet_crm")
}

# The list of {taxon, resource, direction} rows a caveat carries, as a data frame.
#' @noRd
caveats_rows <- function(rows) {
    if (!length(rows)) {
        return(data.frame(taxon = character(0), resource = character(0), direction = character(0),
                          stringsAsFactors = FALSE))
    }
    data.frame(taxon = vapply(rows, function(p) chr(p$taxon), character(1)),
               resource = vapply(rows, function(p) chr(p$resource), character(1)),
               direction = vapply(rows, function(p) chr(p$direction), character(1)),
               stringsAsFactors = FALSE)
}

# The per-value cautions, as a data frame (taxon, resource, cautions, notes).
#' @noRd
cautions_rows <- function(rows) {
    if (!length(rows)) {
        return(data.frame(taxon = character(0), resource = character(0), cautions = character(0),
                          notes = character(0), stringsAsFactors = FALSE))
    }
    data.frame(taxon = vapply(rows, function(p) chr(p$taxon), character(1)),
               resource = vapply(rows, function(p) chr(p$resource), character(1)),
               cautions = vapply(rows, function(p) paste(chr_vector(p$cautions), collapse = " "), character(1)),
               notes = vapply(rows, function(p) paste(chr_vector(p$notes), collapse = "; "), character(1)),
               stringsAsFactors = FALSE)
}

#' @noRd
count_evidence <- function(x, word) {
    sum(x$evidence_consumed == word, na.rm = TRUE) + sum(x$evidence_produced == word, na.rm = TRUE)
}

#' Print CRM parameters and their caveats
#'
#' The caveats are shown every time, because a cell that was never measured, or a link whose size is
#' unknown, has to reach whoever simulates with it.
#'
#' @param x CRM parameters from [foodnet_listen()] or [foodnet_crm()].
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.foodnet_crm <- function(x, ...) {
    cat(sprintf("CRM parameters from %s %s, derived %s from %s\n",
                x$tool, x$tool_version, x$derived_at, x$source_db))
    cat(sprintf("  %d taxa, %d resources; %s phase; values in %s\n", length(x$taxa), length(x$resources),
                x$phase, x$values))
    cat(sprintf("  growth rates: %d of %d taxa (%s)\n", sum(!is.na(x$growth_rates)), length(x$taxa),
                x$growth_rate_unit))
    if (length(x$caveats$media)) {
        cat("  values from: ", paste(x$caveats$media, collapse = " / "), "\n", sep = "")
    }
    if (x$caveats$incomplete) {
        cat("  INCOMPLETE: records could not be read from mGrowthDB, so these numbers may lack data:\n")
        cat(paste0("   * ", x$caveats$errors, "\n"), sep = "")
    }
    cat("  Read before you simulate:\n")
    if (length(x$caveats$warnings)) {
        cat(paste0("   * ", x$caveats$warnings, "\n"), sep = "")
    }
    flagged <- x$caveats$cautions
    if (NROW(flagged)) {
        moving <- flagged[grepl("still_changing", flagged$cautions), , drop = FALSE]
        if (nrow(moving)) {
            named <- paste(moving$taxon, moving$resource)
            cat(sprintf("   * %d value(s) miss use that went on after growth slowed (still_changing): %s%s\n",
                        nrow(moving), paste(utils::head(named, 6), collapse = ", "),
                        if (length(named) > 6) paste0(" and ", length(named) - 6, " more") else ""))
        }
        tiers <- x$caveats$caution_tiers
        if (length(tiers)) {
            # the most serious tier of each value's cautions (Karoline, 2026-10-07: rank them)
            tier_of <- function(words) {
                found <- tiers[intersect(strsplit(words, " ", fixed = TRUE)[[1]], names(tiers))]
                if (length(found)) min(found) else NA
            }
            tier <- vapply(flagged$cautions, tier_of, numeric(1))
            first <- flagged[!is.na(tier) & tier == 1, , drop = FALSE]
            if (nrow(first)) {
                named <- paste(first$taxon, first$resource)
                cat(sprintf("   * %d value(s) carry a caution that changes what they mean (tier 1), read them first: %s%s\n",
                            nrow(first), paste(utils::head(named, 6), collapse = ", "),
                            if (length(named) > 6) paste0(" and ", length(named) - 6, " more") else ""))
            }
            cat(sprintf("   * %d more carry only cautions of certainty (tier 2); x$caveats$cautions lists every caution.\n",
                        sum(!is.na(tier) & tier == 2)))
        } else {
            cat(sprintf("   * %d value(s) carry cautions: x$caveats$cautions lists them.\n", nrow(flagged)))
        }
    }
    phase <- x$phase_growth_rates
    slow <- if (is.null(phase)) character(0) else
        names(phase)[!is.na(phase) & !is.na(x$growth_rates) & phase > 1.1 * x$growth_rates]
    if (length(slow)) {
        cat(sprintf("   * growth rate below the phase's own mean rate for %s: as_miasim(growth = \"phase_floor\")\n",
                    paste(slow, collapse = ", ")))
    }
    single <- count_evidence(x, "single_replicate")
    if (single) {
        cat(sprintf("   * %d cell(s) rest on one replicate (single_replicate).\n", single))
    }
    not_assayed <- count_evidence(x, "not_assayed")
    if (not_assayed) {
        cat(sprintf("   * %d cell(s) were never assayed (NA): no evidence either way.\n", not_assayed))
    }
    unsure <- count_evidence(x, "inconclusive")
    if (unsure) {
        cat(sprintf("   * %d cell(s) are inconclusive (NA): measured, but not pinned down beyond the limit or inside it.\n",
                    unsure))
    }
    no_phase <- count_evidence(x, "no_phase")
    if (no_phase) {
        cat(sprintf("   * %d cell(s) have no phase (NA): measured, but the culture gave no growth phase.\n",
                    no_phase))
    }
    if (x$caveats$mixed_media) {
        cat("   * values pool every medium: the starting concentrations describe no real medium.\n")
    }
    if (length(x$caveats$whole_run)) {
        cat(sprintf("   * no end of growth was found for %s: their values span the whole run, stationary uptake included.\n",
                    paste(x$caveats$whole_run, collapse = ", ")))
    }
    if (x$caveats$stationary_phase) {
        cat("   * these are amounts from after the end of exponential growth, where cells may still grow, stop or\n")
        cat("     die; a CRM reads every uptake as growth, so building one from them needs allow = \"stationary_phase\".\n")
    }
    if (length(x$caveats$biomass_falls)) {
        cat(sprintf("   * biomass fell by more than half over the phase for %s: lysis and death release and take up\n",
                    paste(x$caveats$biomass_falls, collapse = ", ")))
        cat("     compounds, so their amounts need not be the living cells'.\n")
    }
    presence <- x$caveats$presence_only
    if (nrow(presence)) {
        shown <- utils::head(presence, 6)
        cat(sprintf("   * %d link(s) were seen only in another medium, so their size is unknown (NA):\n",
                    nrow(presence)))
        cat("       ", paste0(shown$taxon, " ", shown$direction, " ", shown$resource, collapse = ", "),
            if (nrow(presence) > 6) paste0(" and ", nrow(presence) - 6, " more") else "", "\n", sep = "")
    }
    if (length(x$caveats$conflicts)) {
        cat(sprintf("   * %d value(s) have experiments that disagree (NA): crm_readme(x) names them.\n",
                    length(x$caveats$conflicts)))
    }
    shrunk <- names(x$biomass_change)[!is.na(x$biomass_change) & x$biomass_change <= 0]
    if (length(shrunk) && !x$caveats$stationary_phase) {
        cat(sprintf("   * no biomass gain for %s: a miaSim simulation cannot scale it; leave it out (crm_subset)\n",
                    paste(shrunk, collapse = ", ")))
    }
    if (length(x$caveats$without_a_rate)) {
        cat(sprintf("   * %d taxon(s) have no growth rate: %s\n", length(x$caveats$without_a_rate),
                    paste(x$caveats$without_a_rate, collapse = ", ")))
        cat("       a simulation needs one from elsewhere; as_miasim() stops until you give it.\n")
    }
    if (length(x$other_phases)) {
        other <- names(x$other_phases)[1]
        cat(sprintf("   * the %s phase came too (%s): crm_phase(x, \"%s\").\n", other,
                    if (other == "stationary") "what changed after the end of exponential growth" else "the growth phase",
                    other))
    } else if (x$caveats$searched_both_phases && x$phase == "exponential") {
        cat("   * the search asked for both phases; a CRM describes growth, so these are the exponential phase.\n")
    }
    if (all(is.na(x$biomass_change))) {
        cat("   * no biomass changes came with these parameters, so crm_efficiency() can build only its\n")
        cat("     \"shares\" and \"none\" scales (update foodnet for miaSim's yields).\n")
    }
    cat("   * the cells are measured amounts (net changes), not efficiencies: crm_efficiency(x) turns them\n")
    cat("     into miaSim's E (yields per mM taken up, by-products per unit of growth); crm_backcheck()\n")
    cat("     simulates each taxon alone and says how close it comes to its own monoculture.\n")
    cat("  crm_consumed(x), crm_produced(x), crm_rates(x), crm_resources(x); x$interval_start and x$interval_end\n")
    cat("  hold the hours each value was measured over, x$consumed_lower, x$consumed_upper, x$produced_lower and\n")
    cat("  x$produced_upper its bounds in mM; crm_readme(x) for the full text.\n")
    invisible(x)
}

#' Summarize CRM parameters
#'
#' @param object CRM parameters from [foodnet_listen()] or [foodnet_crm()].
#' @param ... Ignored.
#' @return A list with the counts and the caveats, invisibly; printing it shows the same as [print()].
#' @export
summary.foodnet_crm <- function(object, ...) {
    print(object)
    invisible(list(taxa = length(object$taxa), resources = length(object$resources),
                   with_a_rate = sum(!is.na(object$growth_rates)),
                   not_assayed = count_evidence(object, "not_assayed"),
                   presence_only = nrow(object$caveats$presence_only),
                   conflicts = length(object$caveats$conflicts),
                   without_a_rate = object$caveats$without_a_rate))
}
