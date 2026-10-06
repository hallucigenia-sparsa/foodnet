test_that("the efficiency matrix divides each row by the taxon's uptake", {
    crm <- foodnet:::as_foodnet_crm(example_payload())
    expect_warning(E <- crm_efficiency(crm), "B produced acetate")
    # A: consumed 8 glucose, produced 4 acetate; total uptake 8, so glucose 8/8 = 1, acetate -4/8 = -0.5
    expect_equal(unname(E["A", ]), c(1, -0.5, 0))
    # B: consumed 2 glucose and 6 acetate (total 8) and produced 3 butyrate: 0.25, 0.75, -3/8 = -0.375
    expect_equal(unname(E["B", ]), c(0.25, 0.75, -0.375))
    expect_equal(sum(E["B", E["B", ] > 0]), 1)
})

test_that("without normalizing, E is consumed minus produced in mM", {
    crm <- foodnet:::as_foodnet_crm(example_payload())
    E <- suppressWarnings(crm_efficiency(crm, normalize = "none"))
    expect_equal(unname(E["A", ]), c(8, -4, 0))
})

test_that("NA cells can be refused", {
    crm <- foodnet:::as_foodnet_crm(example_payload())
    expect_error(crm_efficiency(crm, na = "stop"), "NA")
})

test_that("as_miasim stops without a growth rate and builds miaSim's arguments with one", {
    crm <- foodnet:::as_foodnet_crm(example_payload())
    expect_error(as_miasim(crm), "no growth rate for B")
    args <- suppressWarnings(as_miasim(crm, missing_rate = 0.3))
    expect_equal(args$n_species, 2)
    expect_equal(args$n_resources, 3)
    expect_equal(args$names_resources, c("glucose", "acetate", "butyrate"))
    expect_equal(args$growth_rates, c(0.4, 0.3))
    expect_equal(args$resources, c(10, 2, 0))
    expect_equal(dim(args$E), c(2, 3))
})

test_that("the parameters can be written as files", {
    crm <- foodnet:::as_foodnet_crm(example_payload())
    dir <- tempfile()
    on.exit(unlink(dir, recursive = TRUE))
    paths <- crm_write(crm, dir)
    expect_true(all(file.exists(paths)))
    expect_equal(read.csv(paths[3])$growth_rate, c(0.4, NA))
})

test_that("a link that arrives with an amount from another medium is not warned about", {
    payload <- example_payload()
    payload$produced[[2]][[2]] <- 1.5                  # B's acetate, seen elsewhere, now with its value
    crm <- foodnet:::as_foodnet_crm(payload)
    expect_no_warning(E <- crm_efficiency(crm))
    # B: consumed 2 + 6 = 8, produced acetate 1.5 and butyrate 3: acetate 0.75 - 1.5/8 = 0.5625
    expect_equal(unname(E["B", "acetate"]), 0.5625)
})
