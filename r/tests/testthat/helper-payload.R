# A small payload, the shape foodnet.matrix.crm_payload writes, with numbers chosen to be checked by hand.
# Taxon A consumed 8 mM glucose and produced 4 mM acetate; B consumed 2 mM glucose and 6 mM acetate and
# produced 3 mM butyrate; A's butyrate was never assayed; B producing acetate was seen only in another
# medium. B has no growth rate.
example_payload <- function() {
    list(format = "foodnet.crm/v0", tool = "foodnet", tool_version = "0.1.0",
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
         caveats = list(presence_only = list(list(taxon = "B", resource = "acetate", direction = "produced",
                                                  media = list("mMCB"))),
                        conflicts = list(), duplicates = list(), without_a_rate = list("B"),
                        media = list("Wilkins-Chalgren Anaerobe Broth"), value_rule = "majority",
                        searched_both_phases = FALSE),
         readme = "README text\n", studies = list("SMGDB00000001"), settings = list())
}
