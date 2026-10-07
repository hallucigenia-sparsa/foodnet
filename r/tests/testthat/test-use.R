# The example payload (helper-payload.R): A consumed 8 mM glucose, produced 4 mM acetate, grew by 2 at 0.4/h;
# B consumed 2 glucose and 6 acetate, produced 3 butyrate, grew by 1, and has no rate (0.5 where one is given).

with_rate <- function() {
    payload <- example_payload()
    payload$growth_rates <- list(0.4, 0.5)
    foodnet:::as_foodnet_crm(payload)
}

test_that("E follows miaSim: shares of uptake on consumed resources, by-products scaled to the uptake", {
    crm <- with_rate()
    E <- suppressWarnings(crm_efficiency(crm, na = "zero"))
    # A consumed glucose only (8; S = 64) and made 4 acetate: 1, and -4 * 8 / 64
    expect_equal(unname(E["A", ]), c(1, -0.5, 0))
    # B consumed glucose 2 and acetate 6 (C = 8, S = 40) and made 3 butyrate: 0.25, 0.75, -3 * 8 / 40
    expect_equal(unname(E["B", ]), c(0.25, 0.75, -0.6))
    # each taxon's unit: dx * C / (mu * S); A 2 * 8 / (0.4 * 64), B 1 * 8 / (0.5 * 40)
    expect_equal(unname(crm_scale(crm)), c(0.625, 0.4))
})

test_that("a taxon can be left out, and unscaling checks the shape", {
    crm <- crm_subset(with_rate(), taxa = "A")
    expect_equal(crm$taxa, "A")
    expect_equal(dim(crm$consumed), c(1, 3))
    expect_equal(nrow(crm$caveats$presence_only), 0)
    expect_error(crm_unscale(with_rate(), matrix(1, 3, 2)), "one row per taxon")
})

test_that("miaSim's units need a rate and a biomass change for every taxon", {
    crm <- foodnet:::as_foodnet_crm(example_payload())
    expect_error(suppressWarnings(crm_efficiency(crm, na = "zero")), "B \\(no growth rate\\) lacks them")
})

test_that("the shares scale is 0.1.0's matrix, and none is consumed minus produced", {
    crm <- foodnet:::as_foodnet_crm(example_payload())
    E <- suppressWarnings(crm_efficiency(crm, scale = "shares", na = "zero"))
    expect_equal(unname(E["A", ]), c(1, -0.5, 0))
    expect_equal(unname(E["B", ]), c(0.25, 0.75, -0.375))
    E <- suppressWarnings(crm_efficiency(crm, scale = "none", na = "zero"))
    expect_equal(unname(E["A", ]), c(8, -4, 0))
    # miaSim's signs: a resource taken up is positive and a by-product negative, the opposite of foodnet's
    # signed matrix (A consumed 8 mM glucose: -8 there, +8 here; it produced 4 mM acetate: +4 there, -4 here)
    E <- suppressWarnings(crm_efficiency(crm_subset(crm, taxa = "A"), na = "zero"))
    expect_true(E["A", "glucose"] > 0 && E["A", "acetate"] < 0)
})

test_that("every NA set to 0 is counted by its kind", {
    crm <- foodnet:::as_foodnet_crm(example_payload())
    expect_warning(crm_efficiency(crm, scale = "none", na = "zero"), "2 not_assayed, 1 presence_only")
    expect_error(crm_efficiency(crm), "NA")          # the default: NA is never zero unless you say so
})

test_that("booleans are no amounts", {
    payload <- example_payload()
    payload$values <- "booleans"
    expect_error(crm_efficiency(foodnet:::as_foodnet_crm(payload)), "booleans")
})

test_that("as_miasim needs starting abundances and Monod constants, never drawn at random", {
    crm <- with_rate()
    expect_error(as_miasim(crm), "x0")
    expect_error(as_miasim(crm, x0 = crm$biomass_start), "monod_constant")
    args <- suppressWarnings(as_miasim(crm, x0 = crm$biomass_start, monod_constant = 1, na = "zero"))
    expect_equal(args$n_species, 2)
    expect_equal(args$names_resources, c("glucose", "acetate", "butyrate"))
    expect_equal(args$growth_rates, c(0.4, 0.5))
    expect_equal(args$x0, c(0.5 / 0.625, 0.2 / 0.4))         # into each taxon's unit (crm_scale)
    expect_equal(crm_unscale(crm, args$x0), c(A = 0.5, B = 0.2))
    expect_equal(dim(args$monod_constant), c(2, 3))
    expect_equal(args$migration_p, 0)            # no random immigration
    expect_equal(args$resources, c(10, 2, 0))
})

test_that("as_miasim stops without a growth rate", {
    crm <- foodnet:::as_foodnet_crm(example_payload())
    expect_error(as_miasim(crm, x0 = c(1, 1), monod_constant = 1, E = matrix(0, 2, 3)), "no growth rate for B")
})

test_that("pooled media and the stationary phase are refused unless allowed", {
    payload <- example_payload()
    payload$growth_rates <- list(0.4, 0.5)
    payload$caveats$mixed_media <- TRUE
    crm <- foodnet:::as_foodnet_crm(payload)
    expect_error(as_miasim(crm, x0 = c(1, 1), monod_constant = 1), "every medium")
    args <- suppressWarnings(as_miasim(crm, x0 = c(1, 1), monod_constant = 1, allow = "mixed_media", na = "zero"))
    expect_type(args, "list")
    payload$caveats$mixed_media <- FALSE
    payload$caveats$stationary_phase <- TRUE
    expect_error(as_miasim(foodnet:::as_foodnet_crm(payload), x0 = c(1, 1), monod_constant = 1), "stationary")
})

test_that("a taxon simulated alone gains its biomass and makes its by-products when its uptake matches", {
    skip_if_not_installed("miaSim")
    payload <- example_payload()
    payload$growth_rates <- list(0.4, 0.5)
    # A alone, on glucose only, with 8 of 10 mM taken up: biomass and acetate follow the uptake
    payload$taxa <- list("A")
    payload$consumed <- list(list(8, 0, 0))
    payload$produced <- list(list(0, 4, 0))
    payload$evidence_consumed <- list(list("measured", "below_limit", "below_limit"))
    payload$evidence_produced <- list(list("below_limit", "measured", "below_limit"))
    payload[c("growth_rates", "biomass_change", "biomass_start", "biomass_unit", "phase_hours")] <-
        list(list(0.4), list(2), list(0.5), list("OD"), list(12))
    payload$caveats$presence_only <- list()
    payload$caveats$without_a_rate <- list()
    b <- crm_backcheck(foodnet:::as_foodnet_crm(payload), monod_constant = 5)
    taken <- b$simulated[b$what == "glucose"]
    # whatever the uptake, biomass and acetate are in the measured proportion to it
    expect_equal(b$simulated[b$what == "biomass"] / taken, 2 / 8, tolerance = 0.02)
    expect_equal(b$simulated[b$what == "acetate"] / taken, 4 / 8, tolerance = 0.02)
    # and it grows at its measured rate: from 0.5 to 2.5 at 0.4/h takes about log(5) / 0.4 = 4 h when
    # glucose stays saturating (K = 5 mM slows it), never minutes
    hours <- b$simulated[b$what == "hours to grow"]
    expect_gt(hours, 3)
    expect_lt(hours, 12)
})

test_that("a growth rate below the phase's own mean rate is named, and can be floored", {
    payload <- example_payload()
    payload$growth_rates <- list(0.4, 0.5)
    payload$phase_growth_rates <- list(0.8, 0.5)          # A's curves grew faster over the phase than its rate
    crm <- foodnet:::as_foodnet_crm(payload)
    expect_warning(crm_scale(crm), "below the phase's own mean rate for A")
    expect_equal(unname(crm_scale(crm, "phase_floor")), c(2 * 8 / (0.8 * 64), 1 * 8 / (0.5 * 40)))
    args <- suppressWarnings(as_miasim(crm, x0 = crm$biomass_start, monod_constant = 1, na = "zero",
                                       growth = "phase_floor"))
    expect_equal(args$growth_rates, c(0.8, 0.5))
    # a phase rate only a little above the fitted one changes nothing anywhere (the scale's own rule)
    payload$phase_growth_rates <- list(0.42, 0.5)
    crm <- foodnet:::as_foodnet_crm(payload)
    args <- suppressWarnings(as_miasim(crm, x0 = crm$biomass_start, monod_constant = 1, na = "zero",
                                       growth = "phase_floor"))
    expect_equal(args$growth_rates, c(0.4, 0.5))
    # parameters saved by 0.1.0 have no phase rates: the floor changes nothing and loses nothing
    crm$phase_growth_rates <- NULL
    args <- suppressWarnings(as_miasim(crm, x0 = crm$biomass_start, monod_constant = 1, na = "zero",
                                       growth = "phase_floor"))
    expect_equal(args$growth_rates, c(0.4, 0.5))
})

test_that("subsetting keeps the warnings about kept taxa only, once marked, and never an empty one", {
    payload <- example_payload()
    payload$taxa <- list("A", "A b")
    payload$caveats$warnings <- list("A has more data elsewhere", "A b grew slowly", "3 values are inconclusive")
    crm <- foodnet:::as_foodnet_crm(payload)
    kept <- crm_subset(crm, taxa = "A b")
    expect_equal(kept$caveats$warnings, c("(for the whole search) A b grew slowly",
                                          "(for the whole search) 3 values are inconclusive"))
    again <- crm_subset(kept, taxa = "A b")
    expect_equal(again$caveats$warnings, kept$caveats$warnings)
    payload$caveats$warnings <- list("A has more data elsewhere")
    expect_length(crm_subset(foodnet:::as_foodnet_crm(payload), taxa = "A b")$caveats$warnings, 0)
    # the common case: the removed taxon shares no prefix with a kept one
    payload$taxa <- list("Bacteroides fragilis", "Roseburia intestinalis")
    payload$caveats$warnings <- list("no growth curve for Roseburia intestinalis", "2 values are inconclusive")
    kept <- crm_subset(foodnet:::as_foodnet_crm(payload), taxa = "Bacteroides fragilis")
    expect_equal(kept$caveats$warnings, "(for the whole search) 2 values are inconclusive")
    # a removed name inside a kept one, not at its start
    payload$taxa <- list("Clostridium", "[Clostridium] scindens")
    payload$caveats$warnings <- list("no rate for [Clostridium] scindens", "Clostridium has more data elsewhere")
    kept <- crm_subset(foodnet:::as_foodnet_crm(payload), taxa = "[Clostridium] scindens")
    expect_equal(kept$caveats$warnings, "(for the whole search) no rate for [Clostridium] scindens")
})

test_that("a taxon a simulation cannot use is named with what it lacks, and how to leave it out", {
    crm <- foodnet:::as_foodnet_crm(example_payload())                 # B has no growth rate
    expect_error(crm_scale(crm), "and B \\(no growth rate\\) lacks them")
    expect_error(crm_scale(crm), "crm <- crm_subset\\(crm, taxa = setdiff\\(crm\\$taxa, c\\(\"B\"\\)\\)\\)")
})

test_that("passing E where the parameters belong says so", {
    crm <- foodnet:::as_foodnet_crm(example_payload())
    E <- suppressWarnings(crm_efficiency(crm, scale = "shares", na = "zero"))
    expect_error(crm_backcheck(E, monod_constant = 1), "takes the CRM parameters themselves")
    expect_error(crm_efficiency(E), "not a matrix such as crm_efficiency\\(\\)'s E")
})

test_that("a novice is told what to type, in the names the help uses", {
    crm <- foodnet:::as_foodnet_crm(example_payload())               # B has no growth rate
    # the NA stop points at the argument of the function called, not at crm_efficiency()
    expect_error(as_miasim(crm_subset(crm, taxa = "A"), x0 = 1, monod_constant = 1), "pass na = \"zero\" to the function you called")
    # a taxon a simulation cannot use is named before x0 is judged, with a line that runs as typed
    expect_error(as_miasim(crm, x0 = crm$biomass_start, monod_constant = 1, na = "zero"),
                 "crm <- crm_subset\\(crm, taxa = setdiff\\(crm\\$taxa, c\\(\"B\"\\)\\)\\)")
    out <- paste(capture.output(print(crm)), collapse = "\n")
    expect_match(out, "a miaSim simulation cannot use B \\(no growth rate\\)")
})
