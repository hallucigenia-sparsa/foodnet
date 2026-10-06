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
#' taxa by resources, positive where a taxon consumes a resource (its conversion to biomass) and negative
#' where it produces one (a by-product). This builds one from the measured amounts.
#'
#' With `normalize = "consumed"` (the default), each taxon's row is divided by the total it consumed: the
#' consumed entries become the share of its uptake each resource makes up (they sum to 1), and the produced
#' entries the amount of by-product per unit taken up. With `"none"`, `E` is consumed minus produced, in mM.
#' Either is a modeling choice, and a model may need its own scale; that is why the measured amounts are
#' kept apart from this.
#'
#' `NA` cells have no size. With `na = "zero"` they count as 0, and a warning names the links seen only in
#' another medium, which are real links set to 0 here; `"stop"` refuses until you fill them yourself.
#'
#' @param x CRM parameters from [foodnet_listen()] or [foodnet_crm()].
#' @param normalize `"consumed"` or `"none"`, as above.
#' @param na `"zero"` or `"stop"`, as above.
#' @return A numeric matrix, taxa by resources.
#' @examples
#' \dontrun{
#' E <- crm_efficiency(crm)
#' }
#' @export
crm_efficiency <- function(x, normalize = c("consumed", "none"), na = c("zero", "stop")) {
    stopifnot(inherits(x, "foodnet_crm"))
    normalize <- match.arg(normalize)
    na <- match.arg(na)
    consumed <- x$consumed
    produced <- x$produced
    if (anyNA(consumed) || anyNA(produced)) {
        if (na == "stop") {
            stop_foodnet(sum(is.na(consumed)) + sum(is.na(produced)), " cell(s) are NA (not assayed, or seen ",
                         "only in another medium). Fill them, or use crm_efficiency(x, na = \"zero\").")
        }
        presence <- x$caveats$presence_only
        # only the links that reach R without an amount count as 0; with foodnet's "value from the other
        # medium" they arrive with one
        if (nrow(presence)) {
            missing <- mapply(function(taxon, resource, direction) {
                m <- if (direction == "consumed") x$consumed else x$produced
                is.na(m[taxon, resource])
            }, presence$taxon, presence$resource, presence$direction)
            presence <- presence[missing, , drop = FALSE]
        }
        if (nrow(presence)) {
            warning(nrow(presence), " link(s) seen only in another medium count as 0 here, though they are ",
                    "real: ", paste0(utils::head(presence$taxon, 6), " ", utils::head(presence$direction, 6),
                                     " ", utils::head(presence$resource, 6), collapse = ", "),
                    if (nrow(presence) > 6) " and more" else "", call. = FALSE)
        }
        consumed[is.na(consumed)] <- 0
        produced[is.na(produced)] <- 0
    }
    if (normalize == "none") return(consumed - produced)
    total <- rowSums(consumed)
    none <- names(total)[total <= 0]
    if (length(none)) {
        warning("no measured uptake for ", paste(none, collapse = ", "), ": their rows are 0, so a CRM ",
                "gives them nothing to grow on", call. = FALSE)
    }
    scale <- ifelse(total > 0, total, 1)
    (consumed - produced) / scale
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
#' their names, the efficiency matrix `E`, the initial resource concentrations and the growth rates. miaSim
#' is not needed to call this, and no other simulator is assumed: the result is a plain list, so
#' `do.call(miaSim::simulateConsumerResource, c(as_miasim(crm), list(t_end = 48, t_store = 100)))` runs it,
#' and any other model can take the same pieces. (miaSim stores `t_store` time points, 1000 by default, and
#' fails when that exceeds the number of its steps, as it does for a short `t_end`.)
#'
#' A simulation needs a growth rate for every taxon and a starting concentration for every resource, so
#' this stops when one is missing rather than passing a number nobody measured.
#'
#' @param x CRM parameters from [foodnet_listen()] or [foodnet_crm()].
#' @param E An efficiency matrix, by default [crm_efficiency()] of `x`.
#' @param x0 Starting abundances of the taxa, or `NULL` to leave them to the simulator.
#' @param monod_constant A Monod constant matrix (taxa by resources), or `NULL` to leave it to the
#'   simulator. foodnet measures none.
#' @param missing_rate A growth rate for the taxa that have none, or `NULL` to stop when any does.
#' @param missing_resource A starting concentration for the resources that have none, or `NULL` to stop.
#' @return A list with `n_species`, `n_resources`, `names_species`, `names_resources`, `E`, `resources`,
#'   `growth_rates`, and `x0` and `monod_constant` when given.
#' @examples
#' \dontrun{
#' args <- as_miasim(crm)
#' tse <- do.call(miaSim::simulateConsumerResource, c(args, list(t_end = 48, t_store = 100)))
#' }
#' @export
as_miasim <- function(x, E = NULL, x0 = NULL, monod_constant = NULL, missing_rate = NULL,
                      missing_resource = NULL) {
    stopifnot(inherits(x, "foodnet_crm"))
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
    args <- list(n_species = length(x$taxa), n_resources = length(x$resources), names_species = x$taxa,
                 names_resources = x$resources, E = E, resources = unname(resources),
                 growth_rates = unname(rates))
    if (!is.null(x0)) args$x0 <- x0
    if (!is.null(monod_constant)) args$monod_constant <- monod_constant
    args
}
