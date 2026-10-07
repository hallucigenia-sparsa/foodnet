test_that("a payload becomes matrices with names, NA kept as NA", {
    crm <- foodnet:::as_foodnet_crm(example_payload())
    expect_s3_class(crm, "foodnet_crm")
    expect_equal(crm_consumed(crm)["B", "acetate"], 6)
    expect_true(is.na(crm_consumed(crm)["A", "butyrate"]))        # never assayed: NA, not 0
    expect_true(is.na(crm_produced(crm)["B", "acetate"]))         # presence only: NA, not 0
    expect_equal(crm$evidence_produced["B", "acetate"], "presence_only")
    expect_equal(unname(crm_rates(crm)), c(0.4, NA))
    expect_equal(unname(crm_rates(crm, missing = 0.1)), c(0.4, 0.1))
    expect_equal(unname(crm_resources(crm)), c(10, 2, 0))
})

test_that("printing names every caveat", {
    crm <- foodnet:::as_foodnet_crm(example_payload())
    out <- paste(capture.output(print(crm)), collapse = "\n")
    expect_match(out, "2 cell\\(s\\) were never assayed")           # A butyrate, both matrices
    expect_match(out, "1 link\\(s\\) were seen only in another medium")
    expect_match(out, "B produced acetate")
    expect_match(out, "a miaSim simulation cannot use B \\(no growth rate\\)")
})

test_that("an unknown format is read with a warning, a newer major one refused", {
    payload <- example_payload()
    payload$format <- "foodnet.crm/v1-test"
    expect_warning(foodnet:::as_foodnet_crm(payload), "foodnet.crm/v1-test")
    payload$format <- "foodnet.crm/v9"
    expect_error(foodnet:::as_foodnet_crm(payload), "from a newer foodnet")
})

test_that("every value carries the hours it was measured over", {
    x <- foodnet:::as_foodnet_crm(example_payload_both())
    expect_equal(unname(x$interval_start["A", "glucose"]), 0)
    expect_equal(unname(x$interval_end["B", "butyrate"]), 10)
    expect_true(is.na(x$interval_end["A", "butyrate"]))
    # a payload from before 0.2.0 has none: NA, never an error
    old <- foodnet:::as_foodnet_crm(example_payload())
    expect_true(all(is.na(old$interval_start)))
    expect_equal(length(old$other_phases), 0)
})

test_that("with both phases the stationary one can be switched to and back", {
    x <- foodnet:::as_foodnet_crm(example_payload_both())
    expect_equal(x$phase, "exponential")
    expect_output(print(x), "crm_phase\\(x, \"stationary\"\\)")
    s <- crm_phase(x, "stationary")
    expect_equal(s$phase, "stationary")
    expect_true(s$caveats$stationary_phase)
    expect_equal(unname(s$consumed["A", "acetate"]), 1.5)
    expect_equal(unname(s$interval_start["A", "acetate"]), 12)
    expect_equal(unname(s$phase_hours), c(12, 14))
    # the growth rates and the medium are the taxa's, whatever the phase
    expect_identical(s$growth_rates, x$growth_rates)
    expect_identical(s$initial, x$initial)
    back <- crm_phase(s, "exponential")
    expect_identical(back, x)                                  # switching back gives the original
    expect_false(back$caveats$stationary_phase)
    expect_error(crm_phase(foodnet:::as_foodnet_crm(example_payload()), "stationary"), "both phases")
})

test_that("subsetting keeps every phase in step", {
    s <- crm_subset(foodnet:::as_foodnet_crm(example_payload_both()), taxa = "A", resources = c("glucose", "acetate"))
    expect_equal(dim(s$interval_start), c(1, 2))
    expect_equal(dim(s$other_phases$stationary$consumed), c(1, 2))
    expect_equal(unname(crm_phase(s, "stationary")$consumed["A", "acetate"]), 1.5)
})

test_that("each amount carries its bounds, through a change of phase and a subset", {
    x <- foodnet:::as_foodnet_crm(example_payload_both())
    expect_equal(unname(x$consumed_lower["A", "glucose"]), 7.5)
    expect_equal(unname(x$produced_upper["B", "butyrate"]), 3.5)
    expect_true(is.na(x$produced_lower["B", "acetate"]))         # presence only: no number, no bound
    expect_equal(unname(crm_phase(x, "stationary")$consumed_lower["A", "acetate"]), 1.1)
    expect_equal(dim(crm_subset(x, taxa = "B")$consumed_upper), c(1, 3))
    expect_true(all(is.na(foodnet:::as_foodnet_crm(example_payload())$consumed_lower)))
    dir <- tempfile()
    paths <- crm_write(x, dir)
    expect_true(all(file.exists(paths)))
    expect_true(file.exists(file.path(dir, "bounds.csv")) && file.exists(file.path(dir, "intervals.csv")))
})

test_that("stationary amounts are refused by every builder unless allowed", {
    s <- crm_phase(foodnet:::as_foodnet_crm(example_payload_both()), "stationary")
    expect_error(crm_efficiency(s, na = "zero"), "after the end of exponential growth")
    expect_error(crm_efficiency(s, na = "zero"), "allow = \"stationary_phase\"")
    # "shares" and "none" are arithmetic for other models, as in 0.1.0
    expect_true(is.matrix(suppressWarnings(crm_efficiency(s, scale = "none", na = "zero"))))
    expect_error(crm_backcheck(s, monod_constant = 1, na = "zero"), "after the end of exponential growth|miaSim")
    # allowed, it builds (the warning counts the NA cells set to 0)
    expect_true(is.matrix(suppressWarnings(crm_efficiency(s, scale = "shares", na = "zero", allow = "stationary_phase"))))
})

test_that("a switched phase says so when printed and written", {
    p <- example_payload_both()
    p$other_phases$stationary$biomass_falls <- list("A")
    p$other_phases$stationary$biomass_change <- list(-2, 0)
    s <- crm_phase(foodnet:::as_foodnet_crm(p), "stationary")
    out <- paste(capture.output(print(s)), collapse = "\n")
    expect_match(out, "biomass fell by more than half over the phase for A")
    expect_match(out, "the exponential phase came too \\(the growth phase\\)")
    expect_false(grepl("no biomass gain", out))                 # not advice to drop taxa in this phase
    expect_false(grepl("B produced acetate", out))              # the exponential phase's presence-only link
    dir <- tempfile()
    paths <- crm_write(s, dir)
    expect_true(file.exists(file.path(dir, "consumed_stationary.csv")))
    expect_true(file.exists(file.path(dir, "biomass_stationary.csv")))
    expect_true(file.exists(file.path(dir, "growth_rates.csv")))           # not phase files: named as downloaded
    expect_match(readLines(file.path(dir, "README_stationary.txt"))[1], "stationary phase")
    bounds <- utils::read.csv(file.path(dir, "bounds_stationary.csv"))
    expect_true(all(bounds$phase == "stationary"))
    # both phases in one folder: the second write leaves the first's files alone (a review)
    crm_write(crm_phase(s, "exponential"), dir)
    expect_match(readLines(file.path(dir, "README.txt"))[1], "exponential phase")
    expect_true(any(utils::read.csv(file.path(dir, "intervals.csv"))$phase == "exponential"))
    expect_true(file.exists(file.path(dir, "intervals_stationary.csv")))
    expect_equal(bounds$upper_mM[bounds$taxon == "A" & bounds$resource == "acetate" & bounds$direction == "consumed"], 1.9)
})

test_that("second-window compounds keep their window label through a subset", {
    p <- example_payload_both()
    p$resource_phases <- list("exponential", "window", "exponential")       # acetate in the second window
    x <- foodnet:::as_foodnet_crm(p)
    s <- crm_subset(x, resources = c("acetate", "glucose"))
    expect_equal(unname(s$resource_phases), c("window", "exponential"))
    dir <- tempfile()
    crm_write(s, dir)
    rows <- utils::read.csv(file.path(dir, "intervals.csv"))
    expect_equal(unique(rows$phase[rows$resource == "acetate"]), "window")
    expect_equal(unique(rows$phase[rows$resource == "glucose"]), "exponential")
})

test_that("printing carries the page's ranking once, not a list of its own", {
    p <- example_payload()
    p$caveats$caution_tiers <- list(still_changing = 1, pair_decided = 2, growth_rate_boundary = 3)
    p$caveats$warnings <- list("1 change(s) span the wrong or an uneven stretch of time (tier 1; read these first), largest first: A glucose -8.00 mM (exponential: still_changing).")
    p$caveats$cautions <- list(list(taxon = "A", resource = "glucose", cautions = list("growth_rate_boundary", "still_changing"), notes = list()))
    out <- paste(capture.output(print(foodnet:::as_foodnet_crm(p))), collapse = "\n")
    expect_equal(lengths(regmatches(out, gregexpr("tier 1", out))), 1)
    expect_match(out, "x\\$caveats\\$caution_tiers ranks them")
})
