# The resources' chemistry and each taxon's electron balance (foodnet 0.3.0, Karoline 2026-10-08). In the
# example payload (helper-payload.R) glucose has 24 electrons, acetate 8, butyrate 20.

with_chemistry <- function(p = example_payload_both()) {
    p$tool_version <- "0.3.0"
    p$resource_chemistry <- list(
        list(chebi_id = "17234", formula = "C6H12O6", charge = 0, carbon = 6, degree_of_reduction = 24,
             per_cmol = 4, formula_from = "", note = ""),
        list(chebi_id = "30089", formula = "C2H3O2", charge = -1, carbon = 2, degree_of_reduction = 8,
             per_cmol = 4, formula_from = "", note = ""),
        list(chebi_id = "17968", formula = "C4H7O2", charge = -1, carbon = 4, degree_of_reduction = 20,
             per_cmol = 5, formula_from = "", note = ""))
    p$chemistry_source <- list(source = "ChEBI (www.ebi.ac.uk/chebi)", retrieved_at = "2026-10-08")
    p$second_window_resources <- list()
    # A's three cultures, B's two: changes in mM (negative: taken up), NA where a culture measured nothing
    p$culture_changes <- list(
        list(list(culture = "a1", changes = list(-8, 4, NULL)), list(culture = "a2", changes = list(-7.5, 3.6, NULL)),
             list(culture = "a3", changes = list(-8.5, 4.4, NULL))),
        list(list(culture = "b1", changes = list(-2, -6, 3)), list(culture = "b2", changes = list(-2.4, -5.2, 2.5))))
    foodnet:::as_foodnet_crm(p)
}

balance_of <- function(x, taxon) {
    b <- crm_electron_balance(x)
    b[b$taxon == taxon, ]
}

test_that("the chemistry comes as a table, a row per resource", {
    chem <- crm_chemistry(with_chemistry())
    expect_equal(chem$resource, c("glucose", "acetate", "butyrate"))
    expect_equal(chem$degree_of_reduction, c(24, 8, 20))
    expect_equal(chem$charge, c(0, -1, -1))
})

test_that("the electron balance is electrons out over electrons in, by hand", {
    a <- balance_of(with_chemistry(), "A")
    # A: 8 glucose in (192), 4 acetate out (32); its butyrate was never assayed (and the medium holds none)
    expect_equal(a$consumed_electrons_mM, 192)
    expect_equal(a$share, 32 / 192)
    expect_match(a$not_counted, "butyrate: no number")
    expect_equal(a$incomplete, "")
    # the range: the lowest and highest culture's own share (Karoline, 2026-10-08, "Per-culture shares")
    expect_equal(a$cultures, 3)
    expect_equal(a$share_lower, 3.6 * 8 / (7.5 * 24))
    expect_equal(a$share_upper, 4.4 * 8 / (8.5 * 24))
    b <- balance_of(with_chemistry(), "B")
    # B: 2 glucose and 6 acetate in (96), 3 butyrate out (60); its acetate production has no number
    expect_equal(b$share, 60 / 96)
    expect_match(b$not_counted, "acetate: produced: no number")
    expect_equal(c(b$share_lower, b$share_upper), sort(c(60 / 96, 2.5 * 20 / (2.4 * 24 + 5.2 * 8))))
})

test_that("a medium resource without a number is named, and the share stands", {
    p <- example_payload_both()
    p$initial_concentrations <- list(10, 2, 5)
    a <- balance_of(with_chemistry(p), "A")
    expect_equal(a$incomplete, "butyrate")
    expect_equal(a$share, 32 / 192)
})

test_that("a resource without a degree of reduction, or left out, withholds the share", {
    crm <- with_chemistry()
    crm$chemistry$degree_of_reduction[1] <- NA
    a <- balance_of(crm, "A")
    expect_true(is.na(a$share))
    expect_match(a$withheld, "no degree of reduction")
    # a review: crm_subset() without glucose turned the shares into several times their value
    sub <- balance_of(crm_subset(with_chemistry(), resources = c("acetate", "butyrate")), "A")
    expect_true(is.na(sub$share))
    expect_match(sub$withheld, "left out with crm_subset\\(\\) although consumed or produced: glucose")
    # a product left out too (Karoline, 2026-10-08, "Withhold for products too"): A made acetate
    prod <- balance_of(crm_subset(with_chemistry(), resources = c("glucose", "butyrate")), "A")
    expect_true(is.na(prod$share))
    expect_match(prod$withheld, "acetate")
    crm <- with_chemistry()
    crm$chemistry$degree_of_reduction[2] <- NA
    expect_match(balance_of(crm, "A")$withheld, "no degree of reduction")
    # a second-window compound is never counted, so leaving it out changes nothing (a review)
    crm <- with_chemistry()
    crm$second_window_resources <- "acetate"
    left <- crm_subset(crm, resources = c("glucose", "butyrate"))
    expect_equal(balance_of(left, "A")$share, balance_of(crm, "A")$share)
})

test_that("uptake within the detection limits withholds the share", {
    # Karoline, 2026-10-08, "Withhold, say why"; B's stationary phase took up nothing above the limits
    crm <- with_chemistry()
    crm$detection_limits <- c(glucose = 3, acetate = 3, butyrate = 3)
    a <- balance_of(crm, "A")
    # 192 electrons in, against limits of 3 mM x (24 + 8) = 96: a share; at 6 mM, 192 is within them
    expect_equal(a$share, 32 / 192)
    expect_equal(a$cultures_left_out, 0)
    crm$detection_limits[] <- 6
    a <- balance_of(crm, "A")
    expect_true(is.na(a$share))
    expect_match(a$withheld, "within the detection limits")
    # what it is at least: 32 out over 192 in plus 6 x (24 + 8) hidden, and the lowest culture's (a3: 4.4 x 8 out
    # over 8.5 x 24 in plus the same) (Karoline, 2026-10-08)
    expect_equal(a$share_at_least, min(32 / (192 + 192), 3.6 * 8 / (7.5 * 24 + 192), 4.4 * 8 / (8.5 * 24 + 192)))
    expect_true(is.na(a$share_lower))
    # a1 and a2 took up no more than the 192 the limits hide; a3 took up 204
    expect_equal(a$cultures_left_out, 2)
    expect_true(is.na(a$share_upper))
})

test_that("a share outside its cultures' range is withheld", {
    crm <- with_chemistry()
    # A's third culture took up almost nothing and made much: left out of the range, it drives the share
    crm$culture_changes$A[3, ] <- c(-0.1, 9, NA)
    crm$produced["A", "acetate"] <- 5.7
    a <- balance_of(crm, "A")
    expect_equal(a$cultures_left_out, 1)
    expect_true(is.na(a$share))
    expect_match(a$withheld, "outside its cultures' own")
})

test_that("a resource left out that the taxon has no number for stays named", {
    p <- example_payload_both()
    p$initial_concentrations <- list(10, 2, 5)
    left <- crm_subset(with_chemistry(p), resources = c("glucose", "acetate"))
    a <- balance_of(left, "A")
    expect_equal(a$incomplete, "butyrate (left out with crm_subset())")
    expect_equal(a$share, 32 / 192)
    # a subset that leaves nothing out names nothing (a review: an empty name was added)
    expect_equal(balance_of(crm_subset(with_chemistry(p), taxa = "A"), "A")$incomplete, "butyrate")
    # leaving out what made a share undetermined does not bring it back (a review)
    crm <- with_chemistry(p)
    crm$detection_limits <- c(glucose = 4, acetate = 4, butyrate = 4)
    expect_true(is.na(balance_of(crm, "A")$share))
    left <- crm_subset(crm, resources = c("glucose", "acetate"))
    expect_true(is.na(balance_of(left, "A")$share))
    # nor does leaving out a resource measured as 0 (a review): with A's butyrate a measured 0, it is counted
    # and hides 4 x 20; left out, it still does
    crm <- with_chemistry()
    crm$detection_limits <- c(glucose = 4, acetate = 4, butyrate = 4)
    crm$consumed["A", "butyrate"] <- 0
    crm$produced["A", "butyrate"] <- 0
    expect_match(balance_of(crm, "A")$withheld, "within the detection limits")
    zero <- crm_subset(crm, resources = c("glucose", "acetate"))
    expect_match(balance_of(zero, "A")$withheld, "within the detection limits")
    # a resource it never measured, and the medium did not hold, hides nothing when left out (a review)
    crm <- with_chemistry()
    crm$detection_limits <- c(glucose = 4, acetate = 4, butyrate = 4)
    before <- balance_of(crm, "A")
    after <- balance_of(crm_subset(crm, resources = c("glucose", "acetate")), "A")
    expect_equal(after$share, before$share)
    expect_equal(after$withheld, before$withheld)
    # and one held below the limit keeps what it could hide in the floor (a review)
    q <- example_payload_both()
    q$initial_concentrations <- list(10, 2, 0.1)
    low <- with_chemistry(q)
    low$detection_limits <- c(glucose = 6, acetate = 6, butyrate = 6)
    full <- balance_of(low, "A")
    kept <- balance_of(crm_subset(low, resources = c("glucose", "acetate")), "A")
    expect_equal(kept$share_at_least, full$share_at_least)
    # a resource nobody measured in the medium stays named when it alone is left out (a review)
    q$initial_concentrations <- list(10, 2, NULL)
    none <- with_chemistry(q)
    none$detection_limits <- c(glucose = 6, acetate = 6, butyrate = 6)
    expect_match(balance_of(none, "A")$withheld, "taking butyrate")
    expect_match(balance_of(crm_subset(none, resources = c("glucose", "acetate")), "A")$withheld, "taking butyrate")
    # and leaving out a second-window compound keeps the floor's note that it leaves them out (a review)
    none$second_window_resources <- "butyrate"
    expect_match(balance_of(crm_subset(none, resources = c("glucose", "acetate")), "A")$withheld,
                 "leaving out the second time window's compounds")
    expect_match(full$withheld, "at least")
    expect_match(balance_of(left, "A")$withheld, "\\(208 mM\\)")
})

test_that("the medium is judged at the phase's start", {
    p <- example_payload_both()
    p$initial_concentrations <- list(10, 2, 5)
    # A's cultures held no butyrate when the phase began, though the medium did at its first sample
    p$phase_start_mM <- list(list(10, 2, 0), list(10, 2, 5))
    expect_equal(balance_of(with_chemistry(p), "A")$incomplete, "")
})

test_that("a second-window compound is not counted", {
    crm <- with_chemistry()
    crm$second_window_resources <- "acetate"
    a <- balance_of(crm, "A")
    expect_equal(a$share, 0)
    expect_match(a$not_counted, "acetate: measured over the second time window")
})

test_that("the balance follows a subset of taxa and a switch of phase", {
    crm <- crm_subset(with_chemistry(), taxa = "B")
    expect_equal(crm_electron_balance(crm)$share, 60 / 96)
    stationary <- crm_electron_balance(crm_phase(with_chemistry(), "stationary"))
    # after its exponential phase A took up 1.5 mM acetate (12 electrons) and made nothing
    expect_equal(stationary$consumed_electrons_mM[1], 12)
    expect_equal(stationary$share[1], 0)
    back <- crm_phase(crm_phase(with_chemistry(), "stationary"), "exponential")
    expect_equal(crm_electron_balance(back), crm_electron_balance(with_chemistry()))
})

test_that("booleans, no taxa and duplicate names", {
    p <- example_payload_both()
    p$values <- "booleans"
    crm <- with_chemistry(p)
    expect_true(all(is.na(crm_electron_balance(crm)$share)))
    expect_false(any(grepl("crm_electron_balance", capture.output(print(crm)))))
    expect_equal(nrow(crm_electron_balance(crm_subset(with_chemistry(), taxa = character(0)))), 0)
    # two resources of one name (a record without a ChEBI id named as another) still read (a review)
    p <- example_payload_both()
    p$resources <- list("glucose", "acetate", "acetate")
    expect_equal(nrow(crm_chemistry(with_chemistry(p))), 3)
})

test_that("parameters from before 0.3.0 have no chemistry, and say nothing of it", {
    crm <- foodnet:::as_foodnet_crm(example_payload())
    expect_true(all(is.na(crm_chemistry(crm)$degree_of_reduction)))
    expect_true(all(is.na(crm_electron_balance(crm)$share)))
    expect_false(any(grepl("crm_electron_balance", capture.output(print(crm)))))
    expect_true(any(grepl("crm_electron_balance", capture.output(print(with_chemistry())))))
    dir <- tempfile()
    crm_write(crm, dir)
    expect_false(file.exists(file.path(dir, "chemistry.csv")))
})

test_that("the chemistry and the balance are written with foodnet's columns", {
    dir <- tempfile()
    crm_write(with_chemistry(), dir)
    written <- utils::read.csv(file.path(dir, "electron_balance.csv"))
    expect_equal(names(written), c("taxon", "phase", "consumed_electrons_mM", "produced_electrons_mM", "share",
                                   "share_lower", "share_upper", "share_at_least", "cultures", "cultures_left_out",
                                   "withheld", "incomplete", "not_counted"))
    expect_true(file.exists(file.path(dir, "chemistry.csv")))
    crm_write(crm_phase(with_chemistry(), "stationary"), dir)
    expect_true(file.exists(file.path(dir, "electron_balance_stationary.csv")))
})

test_that("R's balance is foodnet's, on parameters foodnet wrote", {
    # crm_synthetic.json is written by foodnet's Python tests (tests/test_chemistry.py checks it is current)
    path <- test_path("crm_synthetic.json")
    crm <- foodnet_crm(path)
    sent <- jsonlite::fromJSON(path, simplifyVector = FALSE)
    for (phase in c("exponential", "stationary")) {
        x <- crm_phase(crm, phase)
        mine <- crm_electron_balance(x)
        theirs <- if (phase == "exponential") sent$electron_balance else sent$other_phases$stationary$electron_balance
        for (i in seq_along(theirs)) {
            for (field in c("consumed_electrons_mM", "produced_electrons_mM", "share", "share_lower",
                            "share_upper", "share_at_least", "cultures", "cultures_left_out")) {
                expected <- theirs[[i]][[field]]
                expect_equal(as.numeric(mine[[field]][i]), if (is.null(expected)) NA_real_ else as.numeric(expected), tolerance = 1e-9,
                             info = paste(phase, i, field))
            }
            expect_equal(mine$incomplete[i], paste(unlist(theirs[[i]]$incomplete), collapse = "; "))
            expect_equal(mine$not_counted[i], paste(unlist(theirs[[i]]$not_counted), collapse = "; "))
            expected <- theirs[[i]]$withheld
            expect_equal(mine$withheld[i], if (is.null(expected)) NA_character_ else expected)
        }
    }
})
