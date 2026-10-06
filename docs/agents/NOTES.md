# Agent notes (shared memory)

The handoff file between agent sessions. Read it at the start of a session and update it in the same pull
request as your change. Record what is not obvious from the code, the git history or the issues; edit in
place; date each entry. Tasks go in GitHub Issues, not here. The repository is public and the house style
gate applies.

## Current state

- 2026-10-04: 0.1.0 built in one session by Karoline's agent from Karoline's specification (her decisions
  are in [docs/METHOD_NOTES.md](../METHOD_NOTES.md), in her words). Public repository
  hallucigenia-sparsa/foodnet created the same day (Karoline's choice of public); CI passed on every job.
  Not released: PyPI trusted publishing and the `pypi` environment are still to set up (RELEASING.md).
- The structure is grow**net** 0.2.0's (crossfeed-bio/crossfeed), copied and adapted: the mGrowthDB client,
  the parallel prefetch, the species list, selection, rates, stats, the page frame, the Cytoscape transport,
  the R transport, the gate and the packaging. The derivation, the model, the matrices, the page content,
  the legend, the help and the R package's content are new.
- Checked against a hand-checked pair of consumption and production matrices (the reference and the
  comparison live outside this repository, since they are unpublished). State at the end of 2026-10-05: with studies 2, 4, 7 and 9 in the second
  box, a 0 to 48 h window and trehalose in the second window, 106 of 108 cells agree as booleans and every
  value of four of the five Wilkins-Chalgren species matches to three decimals. The two remaining cells are
  R. intestinalis acetate and lactate uptake, which the reference measures from each compound's peak to the
  end of the run (a window foodnet deliberately does not offer). The comparison found errors in the
  reference, all corrected there: a duplicate deposit counted twice, two undocumented peak-to-end windows,
  and an F. duncaniae acetate production the data do not support.
- Verified live on 2026-10-04: the page (search, Example, CRM mode, every download), Send to Cytoscape against
  Cytoscape 3.10.3 (node shapes and sizes by kind, arc colors, dashes, widths and phase transparency read back
  from the view), and Send to R through miaSim's `simulateConsumerResource` on the test-case taxa. miaSim
  needs `t_store` below its number of steps for a short `t_end` (documented in `as_miasim`).

- 2026-10-04 (Karoline): the search button reads "Get taxon-metabolite network"; the matrices image comes
  first in the result (`foodnet.figure`, SVG with no dependency; panels stack above `MAX_SIDE_BY_SIDE` px);
  the result table reports the length of the exponential phase. `brand.in_prose` leaves `svg` alone, since
  its name styling is HTML and broke the SVG caption. Tests quote her wording (tests/test_gui.py,
  tests/test_figure.py).

- 2026-10-05 (Karoline): a filled second box is a limit (only matching data, media and ids alike); the
  advanced option `outside_evidence` adds presence from everything else; an empty box considers all data.
  `search.run_query` reads every study of the taxa only when the box holds a medium (it cannot know the
  matching studies before reading), and only the named studies when it holds ids alone. Refined the same
  day: ids are the scope, a medium name the value medium inside it, and with ids alone the majority rule
  runs inside the scope (rule `majority_in_scope`); this reproduces the reference from studies 2, 4, 7, 9.

- 2026-10-05 (Karoline): a second time window for named metabolites (`second_window_*` settings,
  `derive.changes(second=...)`); their cells have phase "window" and `matrix.columns` gives each metabolite
  only the phases it has, labeled by `matrix.interval` when intervals differ. Reproduces every trehalose
  value of the reference (studies 2, 4, 7, 9; 0 to 48 h; trehalose to the last sample).

## Things learned about mGrowthDB while building

- **Search misses culture-level monocultures.** `search.json?strainNcbiIds=` matches per-strain measurement
  contexts only, so study SMGDB00000009 (E. coli LF82 and B. fragilis, culture-level OD and flow cytometry)
  never comes back. `SpeciesIndex.where` (taxon id to studies, from the crawl) is used instead. grow**net**
  may have the same gap for monoculture-only studies.
- **Duplicate deposits.** Study 2 RI_WC equals study 7 ri2 (rounded); study 7 bt2 equals study 2 BT_WC.
  RI_WC also holds one acetate series twice (replicates 1 and 2), so containment is tested one way.
- **Acid and base names.** Studies 4 and 9 use the acids (acetic, butyric, formic, propionic, isobutyric,
  isovaleric, (S)-lactic), the others the bases. `compounds.EQUIVALENT` joins them.
- **Units.** mGrowthDB serves nearly every concentration in mM already (`techniqueUnits`). Study 19 records
  a metabolite in "AUC", which is refused as not a concentration.
- **Two spellings of Wilkins-Chalgren** ("... Broth (WC)" and "... Broth") are one medium key.
- **Study 10** labels experiments "Roseburia intestinalis" and "Lachnospiraceae bacterium 7_1_58FAA" with the
  same strain, and holds F. duncaniae in Db-MM: for Karoline's mGrowthDB list.
- **Zeros that look like missing values** in study 9 (B. fragilis and E. coli, one whole metabolite panel
  at one time point each, plus B. fragilis formate at 4 h) are used as served: Karoline decided against
  treating them as missing (METHOD_NOTES decision 15). Do not add a filter.

## Conventions kept from grow**net**

- A number never computed is missing, never 0 (JSON null, left out of GraphML and the Cytoscape payload,
  NA in the matrices).
- CyREST: styles and layouts are applied with GET; an existing style is updated in place; columns an arc may
  lack are declared after the network is posted.
- The name in running text is styled (`brand.in_prose`); in the READMEs it is written food**net**.
- R: miaSim is suggested, never imported; the listener is base R sockets; `make rd` writes the help pages
  from the `#'` comments (no roxygen2 needed); `make r-check` runs `R CMD check`.
- The R package listens on 8794, grow**net**'s on 8793, so both can run in one R session.
