# Changelog

All notable changes to foodnet. The format follows Keep a Changelog, and the version numbers follow
semantic versioning.

## [Unreleased]

### Added

- `foodnet --version` prints the version.
- The help documents the command line in full: every command and option (read from the command line itself, so
  none can go missing) and which file each output holds, including the tables found in both the matrices zip and
  the CRM zip. `foodnet derive --help` groups its options by what they are for.
- The CRM README, the help and the README say over which interval the growth rates are fitted: each
  replicate's whole growth curve, by the rate method's own interval, not the phase or time window of the values.
- `growth_rates.csv` (and the report, crm.json's `growth_rate_detail` and R's `crm_write()`) gives `fit_start_h`
  and `fit_end_h`: the earliest first and the latest last hour any replicate's growth-rate fit used.
- `growth_rates.csv` and `decline_rates.csv` name the growth curve's technique (technique, unit, strain or
  culture), and the CRM README and the report warn when taxa rest on different techniques: a cell count can
  fall where optical density stays flat.
- Both downloads' README says where the numbers come from: which mGrowthDB records (and which growth curve of
  a replicate), what only ChEBI gives (formula and charge), what foodnet computes from them, and the studies
  read with their uploadedAt.
- CRM mode gives each taxon's decline rate (`decline_rates.csv`, `decline_rates` and `decline_rate_detail` in
  `crm.json`, `crm$decline_rates` in R): the loss of natural-log abundance per hour from the end of growth (the
  first sample at 90% of the linear rise) to the last sample, the mean over the replicates of one technique,
  with the scatter between replicates, the number of experiments, the window, how well a line fits (fit_r2),
  and whether the curves leave the plateau or fall to zero (death) or stay on it. Absent, not 0, where there
  is no usable decline.
- CRM mode gives the Pirt maintenance coefficient as a separate characteristic (`maintenance.csv`,
  `maintenance` in `crm.json`, `crm$maintenance` in R), which nothing else uses: each taxon's uptake of each
  compound per unit of biomass per hour on its growth plateau (from where the decline starts), with its
  scatter, experiments, window, the
  compound's concentration at the plateau's start and how many replicates had it below the detection limit
  then.

### Fixed

- A count of dead cells (a flow cytometry gate) could become a culture's growth curve when the others failed; it
  never does now. The gate of a count (live, total) is part of its technique in every table, and a per-strain
  count that repeats a gated one is named for that gate: study SMGDB00000007's per-strain counts are live counts.
- A duplicate deposit, counted once in the matrices, was counted twice in the growth rates' replicates (and so
  their n and median), in the biomass and in the initial concentrations: it now counts once, and the copy
  with the richer growth curve gives the growth rate.
- Growth rates take the growth curves of the technique and unit most of a taxon's replicates use (they pooled
  OD and cell counts), leave out cultures that did not grow and cultures whose values span the whole run where
  the taxon has others with a phase, and a growth-only monoculture's growth curve is chosen and tested as the
  others are.
- On macOS foodnet checks mGrowthDB's and ChEBI's certificates against the system's trust store, through
  truststore 0.10.5 or later, now its one runtime dependency (macOS only), with a context per thread. With a
  Python uv put on a Mac, which lacks the root mGrowthDB's certificate chains to, `uv tool install foodnet`
  failed at the first search.

## [0.4.0] (2026-10-09)

### Added

- foodnet keeps two things in the user's cache folder: mGrowthDB's species list, for a day, and what it read of
  each study (experiments, replicates, their series), while the study's `uploadedAt` is unchanged and for at
  most 30 days. A repeated search reads only the study records (4 requests for the CRM example instead of 734).
  `foodnet derive --refresh`, or Read everything from mGrowthDB again in Advanced settings, uses neither.
- The data versions of a result (`data` in the network JSON and `crm.json`, and the report) say when the species
  list was read and which studies came from the copy on this machine, since when; a name not found says the
  date of the species list it was looked up in.
- When mGrowthDB answers slowly, the progress text says so: how many requests took over 10 s or were retried.
- The About page lists the releases, and the help has a section on the cautions: what each means, by tier,
  and what to do about it.
- `foodnet derive --taxa` takes names separated by commas or semicolons, so strain names need no quotes
  (`--taxa Escherichia coli LF82, Bacteroides fragilis`).

### Changed

- `foodnet derive --taxa`: once a comma or semicolon appears, the bare words between two of them make one name,
  so `--taxa Blautia Roseburia,Bacteroides` now looks for "Blautia Roseburia" (reported as not found) where 0.3.0
  read three names; a quoted name is always one name.
- A replicate's series are read from its one CSV (`bioreplicate/<id>.csv`), one request per replicate instead
  of one per series (53 instead of 542 for the CRM example); the values are the same.
- The study crawl no longer stops at a gap of absent study ids before the last study it found before; past
  that study, 25 absent ids in a row end it, as before.

### Fixed

- A certificate the Python foodnet runs on cannot verify (some Pythons uv installs on macOS lack the root
  mGrowthDB's chain ends at) is now said so, with what to do, instead of "check your network", and is not
  retried.
- The data versions (species list, studies from foodnet's copy) are in every download's README and in R's
  `print()` and `x$data`; the About page's changelog link names the version.
- A server address that does not accept a connection within 3 s is left for the next one (and tried again with
  the full time if none answers). mGrowthDB's IPv6 address did not answer through the KU Leuven VPN, and every
  new connection waited about 17 s for it, which made a search take several minutes.
- A stationary value pooled with an experiment that reached no end of growth rests on the other experiments
  and names that one in its notes; it carried `no_phase`, which is no caution (its arc failed validation), and
  credited that experiment's study. The legend wraps long caution lines instead of cutting them off,
  and several of them say more plainly what they mean.

## [0.3.0] (2026-10-08)

### Added

- CRM mode reads each resource's formula and charge from ChEBI and gives its carbon atoms and degree of
  reduction (`chemistry.csv`; `resource_chemistry` in `crm.json`), and each taxon's electron balance: the
  electrons in its measured by-products over those in what it consumed, with the lowest and highest culture's
  own share (`electron_balance.csv`; `electron_balance` in `crm.json` and in each of `other_phases`, with each
  culture's changes in `culture_changes`). The report lists both. It is a check and changes no value. It is
  withheld when a resource the taxon consumed or produced has no degree of reduction, when the share lies
  outside its own cultures' range, or when the electrons taken up are within the detection limits of the
  resources counted or missing (then it says what the share is at least, `share_at_least`). It names the resources the medium
  held when the phase began that the taxon has no number for (`phase_start_mM` in `crm.json`, with each
  resource's detection limit in `detection_limits_mM`), and `cultures_left_out` counts the cultures that took
  up no more than those limits hide. `chemistry_source` names ChEBI and the retrieval date. With growth rates
  on, foodnet now also contacts ChEBI's public API; if ChEBI cannot be reached (each lookup is tried twice,
  then ChEBI is left alone for the rest of that search), the parameters still come, with a warning and these
  fields blank.
- R: `crm_chemistry()` and `crm_electron_balance()`, which follow `crm_subset()` (a resource left out that the
  taxon consumed or produced withholds the share) and `crm_phase()`; `crm_write()` writes both files with the
  download's columns, and `print()` points to them.
- `crm.json` names the compounds of the second time window (`second_window_resources`).

### Fixed

- Isovalerate was keyed under CHEBI:50128, which is biflavonoid; it is CHEBI:48942, as mGrowthDB records it.
  Isovaleric acid (CHEBI:28484) now joins isovalerate instead of standing apart. The metabolite's id in the
  network files, matrices and `initial_concentrations.csv` changes from `chebi:50128` to `chebi:48942`; no
  value changes with today's mGrowthDB, where isovalerate (48942) is recorded only in a community.

## [0.2.0] (2026-10-07)

### Upgrading from 0.1.0 (breaking changes)
- **Install the R package of the same version as foodnet**: the R listener now needs a secret the page reads
  from a file, so a 0.1.0 page cannot send to a 0.2.0 package, nor a 0.2.0 page to a 0.1.0 one. The page, the
  help and the READMEs print the install line for the version (`ref = "v0.2.0"`). The R package now refuses
  parameters from a newer major payload format instead of reading them.
- R: `as_miasim(x, x0, monod_constant, E = NULL, ...)` takes the starting abundances and Monod constants
  second and third and requires them (0.1.0: `as_miasim(x, E, x0, monod_constant)`, and miaSim drew them at
  random when missing); `crm_efficiency()` builds miaSim's E in each taxon's own unit (`scale = "miasim"`),
  replacing `normalize =` (`scale = "shares"` is 0.1.0's matrix); `na` defaults to `"stop"`, so NA cells are
  refused until you say how (`na = "zero"`).
- The CRM payload is `foodnet.crm/v1` (0.1.0: v0) and the network format `foodnet.metabolite_network/v1`
  (arcs gained `n_experiments` and the phase `whole_run`, which the v0 schema refuses).
- The first header cell of every matrix CSV names the value medium (`taxon [values from ...]`), so code that
  finds the row names by the header `taxon` must take the first column instead.
- The command line exits 3 when records could not be read (`--allow-incomplete` accepts the result).
- Values are judged differently (below), so cells can move between a number, 0 and NA.

### Added
- CRM parameters say, cell by cell, the hours each value was measured over (`interval_start_h` and
  `interval_end_h` in crm.json, intervals.csv in the download; `x$interval_start` and `x$interval_end` in R):
  the window the change was measured over, the mean over its experiments.
- CRM parameters carry bounds on each amount (`consumed_lower`, `consumed_upper`, `produced_lower`,
  `produced_upper` in crm.json and in R, bounds.csv in the download): the lowest and highest replicate, as
  amounts clipped at 0; a 0 runs from 0 to the detection limit at least.
- With Both, the CRM parameters carry the stationary phase beside the exponential one (`other_phases` in
  crm.json; consumed_stationary.csv, produced_stationary.csv, their evidence and biomass_stationary.csv in the
  download): what changed after the end of exponential growth, with the caveat `biomass_falls` for taxa whose
  biomass fell by more than half. In R, `crm_phase(x, "stationary")` switches to it; `crm_efficiency()`,
  `as_miasim()` and `crm_backcheck()` build from it only with `allow = "stationary_phase"`.
- R: `crm_scale()`, `crm_unscale()`, `crm_subset()`, `crm_backcheck()`, `crm_phase()`, and
  `growth = "phase_floor"`; `crm_write()` writes the download's files.
- Cautions within_scatter, within_evaporation (with the Evaporation setting), pair_decided and start_differs.
- A CRM example button (five gut species in Wilkins-Chalgren, CRM mode on) and a help section that walks through
  a miaSim simulation step by step, which the example runs as it is.
- Cautions are ranked: tier 1, the value spans the wrong or an uneven stretch of time; tier 2, whether there is
  a change, or how large, is less certain; tier 3, how the phase was found. The page's first warning, also in
  the files' README and R's print, names the tier 1 values: first those still changing after the growth rate
  fell, by the share of the change their exponential value misses, then the others by size. still_changing
  tests the whole stretch from the growth-rate boundary to where the 90% rule would have ended growth.
- R: `crm_backcheck()` warns that over a time window neither of its rows judges the Monod constants; the
  help says to choose them on the exponential phase.
- R: each function taking the parameters says so when handed something else (crm_efficiency()'s E, say), and
  a taxon a miaSim simulation cannot use is named with what it lacks and the crm_subset() line that leaves it
  out, by crm_scale(), as_miasim() and print().
- The result opens with its buttons (downloads, Cytoscape, report, CRM parameters), then an index of its
  sections, then the matrices image.

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
- New evidence states in the matrices: inconclusive, no_phase, seen_elsewhere, single_replicate, whole_run
  and not_grown; the image draws them. Every matrix CSV names its value medium (and INCOMPLETE when it is) in its first header cell,
  and the zips hold cautions.csv and the page's warnings in their README.
- Only taxa with phase values vote for the value medium (not whole-run or second-window values, so naming a
  compound for the second window never changes the medium); the page names taxa with more data in another
  medium, a narrow win, and taxa without a growth phase.
- Merging to genus applies the detection limit and makes disagreeing taxa inconclusive.
- Exponential growth ends at the 90% rule or earlier, where the growth rate has fallen below a tenth of its
  maximum over two consecutive intervals (E. coli LF82 now ends at 8 h, not 84 to 108 h); the first growth
  curve that gives a boundary is used; coarse sampling and replicates that end growth far apart are flagged.
- `crm_efficiency()` and `as_miasim()` refuse NA cells unless told `na = "zero"`; each listener on a port has
  its own secret; a request that has not sent its headers within 5 s is dropped.
- A culture without an end of exponential growth gives its change over the whole run, marked whole_run (arcs
  in phase whole_run, a dashed frame in the image), instead of no value, when it grew by the 1.5-fold rule or
  its optical density rose by 0.1, or when its compounds moved as metabolism moves them (growth_unclear); one
  that did neither gives no value (not_grown). The fold rule reads the curve smoothed by a running median.
- R: as_miasim() and crm_backcheck() switch off miaSim's random immigration (migration_p = 0), which it adds
  even when stochastic is FALSE and which swamped growth in foodnet's units.
- Media are told apart by the amounts their descriptions state (in one notation) and by more phrasings.
- CRM parameters (`foodnet.crm/v1`) carry each taxon's biomass change; the R package builds miaSim's
  efficiency matrix from its equations, requires starting abundances and Monod constants, refuses pooled
  media and the stationary phase unless allowed, and adds `crm_backcheck()`.
- Identical replicates give no test; multiple testing is corrected over the tests that make the arcs.
- Duplicate deposits are matched with rounded sample times and need a series that moves.
- Succinic and valeric acid join their bases.

### Fixed
- Hovering over a cell of the matrices image shows one box, not also the browser's own tooltip.
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
