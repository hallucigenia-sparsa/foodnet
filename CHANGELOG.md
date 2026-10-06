# Changelog

All notable changes to foodnet. The format follows Keep a Changelog, and the version numbers follow
semantic versioning.

## [Unreleased]

### Changed
- A cell's value counts each experiment once: the mean and its test are over the experiments' means;
  `n_experiments` is a new arc field.
- A change is judged on a one-sided 90% confidence interval on the mean against the detection limit with
  three or more replicates (no change: the interval inside it), on both replicates with two, and not at all
  on one; experiments must not contradict each other. Otherwise the value is inconclusive (NA, no arc)
  instead of an arc or a 0. Experiments that agree in direction but not in size
  carry the caution amounts_differ; an inconclusive experiment that does not contradict the others is left
  out and named (experiment_left_out). The end of exponential growth by the growth rate is marked
  (growth_rate_boundary), and a compound still changing right after it (still_changing).
  Advanced settings can judge by the mean alone, and give compounds their own detection limits.
- New evidence states in the matrices: inconclusive, no_phase, seen_elsewhere and single_replicate; the image
  draws them. Every matrix CSV names its value medium (and INCOMPLETE when it is) in its first header cell,
  and the zips hold cautions.csv and the page's warnings in their README.
- Only taxa with values vote for the value medium; the page names taxa with more data in another medium, a
  narrow win, and taxa without a growth phase.
- Merging to genus applies the detection limit and makes disagreeing taxa inconclusive.
- Exponential growth ends at the 90% rule or earlier, where the growth rate has fallen below a tenth of its
  maximum over two consecutive intervals (E. coli LF82 now ends at 8 h, not 84 to 108 h); the first growth
  curve that gives a boundary is used; coarse sampling and replicates that end growth far apart are flagged.
- `crm_efficiency()` and `as_miasim()` refuse NA cells unless told `na = "zero"`; each listener on a port has
  its own secret; a request that has not sent its headers within 5 s is dropped.
- Media are told apart by the amounts their descriptions state (in one notation) and by more phrasings.
- CRM parameters (`foodnet.crm/v1`) carry each taxon's biomass change; the R package builds miaSim's
  efficiency matrix from its equations, requires starting abundances and Monod constants, refuses pooled
  media and the stationary phase unless allowed, and adds `crm_backcheck()`.
- Identical replicates give no test; multiple testing is corrected over the tests that make the arcs.
- Duplicate deposits are matched with rounded sample times and need a series that moves.
- Succinic and valeric acid join their bases.

### Fixed
- The R listener accepted parameters from any machine on the network and from web pages; it now needs a
  one-time secret only this user can read, and a stray request can no longer stop it.
- A result built on records that could not be read gave no sign in its files and exit status 0; it now says
  `incomplete` and the command line exits 3 (`--allow-incomplete` accepts it).
- A species list built with read failures dropped taxa silently and was kept for an hour.
- A download naming a search that was gone served another search's result.
- Study links accepted any scheme; CSV text cells could run as spreadsheet formulas.

## [0.1.0] (2026-10-06)

### Added
- Bipartite taxon and metabolite networks from mGrowthDB batch monocultures: produced arcs from a taxon to a
  metabolite, consumed arcs from a metabolite to a taxon, with the amount in mM as the width.
- Growth phases: Exponential phase (the default), Stationary phase or Both, with the end of exponential
  growth estimated from each replicate's growth curve; a time window in Advanced settings replaces them.
- Values from the medium holding data for the most taxa, or from the media, experiments or studies in the
  second box; other media give presence only. Advanced: ignore media differences, report booleans, the
  detection limit (0.2 mM), merge arcs across studies, merge to genus.
- Duplicate deposits counted once; pooled experiments that disagree reported as conflicts; acid and base
  forms of a compound joined.
- Outputs: JSON, GraphML, a taxa by metabolites matrix, the consumed and produced matrices with their
  evidence, Send to Cytoscape, and the report.
- CRM mode: growth rates and initial medium concentrations, downloadable or sent to R, and an R companion
  package that builds miaSim's consumer-resource arguments.
- The local page, the command line, the Windows program, and grownet's guardrail gate and release checks.
- The consumed and produced matrices as an image (SVG), shown first among the results and downloadable;
  hovering over a cell shows its value and the studies, experiments and medium behind it.
- The length of the exponential phase per arc, in the result table and every output.
- Media are told apart by what their descriptions say was added or taken away and by a recorded atmosphere
  (on by default, switchable); a medium named in the second box gives values from the medium it matches for
  the most taxa, the others it matches presence. Experiments can be excluded by id.
- An advanced setting for matrix cells seen only in another medium: NA (default), TRUE, or the amount
  measured there, drawn on its own background in the image.
- A filled second box limits the search: study or experiment ids set the scope, a medium name the value
  medium, and with ids alone the majority medium within them gives the values and their other media
  presence. The advanced option "Include supporting evidence outside the second box" adds the rest as
  presence. An empty box considers all data.
- A second time window for metabolites named in Advanced settings (an empty end means the last sample), so
  a compound such as trehalose can be measured over the whole run while the rest use the phases or a window.
- Strains are named by their current name in mGrowthDB, as in grownet.
- A second box that matches none of the taxa's monocultures says so and lists the media they were grown in.
- No metabolite is left out by default; "Leave out these metabolites" is an advanced setting.
