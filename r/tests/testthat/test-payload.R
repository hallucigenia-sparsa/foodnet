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
    expect_match(out, "1 taxon\\(s\\) have no growth rate: B")
})

test_that("an unknown format is read with a warning", {
    payload <- example_payload()
    payload$format <- "foodnet.crm/v9"
    expect_warning(foodnet:::as_foodnet_crm(payload), "foodnet.crm/v9")
})
