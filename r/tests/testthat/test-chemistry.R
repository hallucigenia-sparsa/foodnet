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
    foodnet:::as_foodnet_crm(p)
}

test_that("the chemistry comes as a table, a row per resource", {
    chem <- crm_chemistry(with_chemistry())
    expect_equal(chem$resource, c("glucose", "acetate", "butyrate"))
    expect_equal(chem$degree_of_reduction, c(24, 8, 20))
    expect_equal(chem$charge, c(0, -1, -1))
})

test_that("the electron balance is electrons out over electrons in, by hand", {
    b <- crm_electron_balance(with_chemistry())
    # A: 8 glucose in (192), 4 acetate out (32); its butyrate was never assayed
    expect_equal(b["A", "consumed_e_mM"], 192)
    expect_equal(b["A", "share"], 32 / 192)
    expect_match(b["A", "not_counted"], "butyrate: no number")
    # B: 2 glucose and 6 acetate in (96), 3 butyrate out (60); its acetate production has no number
    expect_equal(b["B", "share"], 60 / 96)
    expect_match(b["B", "not_counted"], "acetate: produced: no number")
    # the range: the fewest electrons out over the most in, and the reverse
    expect_equal(b["A", "lower"], 3.6 * 8 / (8.5 * 24 + 0.1 * 8))
    expect_equal(b["A", "upper"], (0.1 * 24 + 4.4 * 8) / (7.5 * 24))
    expect_equal(b["B", "lower"], 2.5 * 20 / (2.4 * 24 + 6.8 * 8 + 0.2 * 20))
    expect_equal(b["B", "upper"], (0.1 * 24 + 3.5 * 20) / (1.6 * 24 + 5.2 * 8))
})

test_that("a second-window compound is not counted", {
    p <- example_payload_both()
    crm <- with_chemistry(p)
    crm$second_window_resources <- "acetate"
    b <- crm_electron_balance(crm)
    expect_equal(b["A", "share"], 0)
    expect_match(b["A", "not_counted"], "acetate: measured over the second time window")
})

test_that("the balance follows a subset and a switch of phase", {
    crm <- crm_subset(with_chemistry(), resources = c("glucose", "butyrate"))
    expect_equal(crm_chemistry(crm)$resource, c("glucose", "butyrate"))
    expect_equal(crm_electron_balance(crm)["B", "share"], 60 / 48)
    stationary <- crm_electron_balance(crm_phase(with_chemistry(), "stationary"))
    # after its exponential phase A took up 1.5 mM acetate (12 electrons) and made nothing
    expect_equal(stationary["A", "consumed_e_mM"], 12)
    expect_equal(stationary["A", "share"], 0)
})

test_that("parameters from before 0.3.0 have no chemistry, and say nothing of it", {
    crm <- foodnet:::as_foodnet_crm(example_payload())
    expect_true(all(is.na(crm_chemistry(crm)$degree_of_reduction)))
    expect_true(all(is.na(crm_electron_balance(crm)$share)))
    expect_false(any(grepl("crm_electron_balance", capture.output(print(crm)))))
    expect_true(any(grepl("crm_electron_balance", capture.output(print(with_chemistry())))))
})

test_that("the chemistry and the balance are written with the other files", {
    dir <- tempfile()
    crm_write(with_chemistry(), dir)
    expect_true(all(file.exists(file.path(dir, c("chemistry.csv", "electron_balance.csv")))))
    crm_write(crm_phase(with_chemistry(), "stationary"), dir)
    expect_true(file.exists(file.path(dir, "electron_balance_stationary.csv")))
})
