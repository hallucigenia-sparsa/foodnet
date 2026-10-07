# A small payload, the shape foodnet.matrix.crm_payload writes, with numbers chosen to be checked by hand.
# Taxon A consumed 8 mM glucose and produced 4 mM acetate; B consumed 2 mM glucose and 6 mM acetate and
# produced 3 mM butyrate; A's butyrate was never assayed; B producing acetate was seen only in another
# medium. B has no growth rate. A grew by 2 (cells/mL, say) from 0.5 over 12 h; B by 1 from 0.2 over 10 h.
example_payload <- function() {
    list(format = "foodnet.crm/v1", tool = "foodnet", tool_version = "0.1.0",
         derived_at = "2026-10-04T12:00:00+02:00", source_db = "test", phase = "exponential",
         values = "mM", detection_limit_mM = 0.2,
         taxa = list("A", "B"), taxon_ids = list("ncbi:1", "ncbi:2"),
         resources = list("glucose", "acetate", "butyrate"),
         resource_ids = list("chebi:17234", "chebi:30089", "chebi:17968"),
         consumed = list(list(8, 0, NULL), list(2, 6, 0)),
         produced = list(list(0, 4, NULL), list(0, NULL, 3)),
         evidence_consumed = list(list("measured", "below_limit", "not_assayed"),
                                  list("measured", "measured", "below_limit")),
         evidence_produced = list(list("below_limit", "measured", "not_assayed"),
                                  list("below_limit", "presence_only", "measured")),
         growth_rates = list(0.4, NULL), growth_rate_unit = "1/h",
         initial_concentrations = list(10, 2, 0), initial_unit = "mM",
         biomass_change = list(2, 1), biomass_start = list(0.5, 0.2), biomass_unit = list("OD", "OD"),
         phase_hours = list(12, 10),
         caveats = list(presence_only = list(list(taxon = "B", resource = "acetate", direction = "produced",
                                                  media = list("mMCB"))),
                        conflicts = list(), duplicates = list(), without_a_rate = list("B"),
                        media = list("Wilkins-Chalgren Anaerobe Broth"), value_rule = "majority",
                        searched_both_phases = FALSE, mixed_media = FALSE, stationary_phase = FALSE,
                        inconclusive = list()),
         readme = "README text\n", studies = list("SMGDB00000001"), settings = list())
}

# The same with both phases (foodnet 0.2.0, Karoline 2026-10-07): every value with the hours it was measured
# over, and the stationary phase beside the exponential one. A's exponential phase runs 0 to 12 h, B's 0 to
# 10 h; afterwards A takes up 1.5 mM acetate to 24 h, without growth.
example_payload_both <- function() {
    p <- example_payload()
    p$tool_version <- "0.2.0"
    p$interval_start_h <- list(list(0, 0, NULL), list(0, 0, 0))
    p$interval_end_h <- list(list(12, 12, NULL), list(10, 10, 10))
    p$consumed_lower <- list(list(7.5, 0, NULL), list(1.6, 5.2, 0))
    p$consumed_upper <- list(list(8.5, 0.1, NULL), list(2.4, 6.8, 0.2))
    p$produced_lower <- list(list(0, 3.6, NULL), list(0, NULL, 2.5))
    p$produced_upper <- list(list(0.1, 4.4, NULL), list(0.1, NULL, 3.5))
    p$phases <- list("exponential", "stationary")
    p$caveats$searched_both_phases <- TRUE
    p$other_phases <- list(stationary = list(
        phase = "stationary",
        consumed = list(list(0, 1.5, NULL), list(0, 0, 0)),
        produced = list(list(0, 0, NULL), list(0, NULL, 0)),
        evidence_consumed = list(list("below_limit", "measured", "not_assayed"),
                                 list("below_limit", "below_limit", "below_limit")),
        evidence_produced = list(list("below_limit", "below_limit", "not_assayed"),
                                 list("below_limit", "presence_only", "below_limit")),
        interval_start_h = list(list(12, 12, NULL), list(10, 10, 10)),
        interval_end_h = list(list(24, 24, NULL), list(24, 24, 24)),
        consumed_lower = list(list(0, 1.1, NULL), list(0, 0, 0)),
        consumed_upper = list(list(0.1, 1.9, NULL), list(0.1, 0.1, 0.1)),
        produced_lower = list(list(0, 0, NULL), list(0, NULL, 0)),
        produced_upper = list(list(0.1, 0.1, NULL), list(0.1, NULL, 0.1)),
        biomass_change = list(0, 0), biomass_start = list(2.5, 1.2), biomass_unit = list("OD", "OD"),
        phase_hours = list(12, 14), phase_growth_rates = list(0, 0),
        cautions = list(), inconclusive = list()))
    p
}
