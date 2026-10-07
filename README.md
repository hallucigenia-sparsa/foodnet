# foodnet: who produces and who consumes which metabolite

[![ci](https://github.com/hallucigenia-sparsa/foodnet/actions/workflows/ci.yml/badge.svg)](https://github.com/hallucigenia-sparsa/foodnet/actions/workflows/ci.yml)
[![license: Apache-2.0](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](LICENSE)
[![python: 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](pyproject.toml)

food**net** builds bipartite taxon and metabolite networks from the batch monocultures in
[mGrowthDB](https://mgrowthdb.gbiomed.kuleuven.be/): which taxon produces which compound, and which consumes
it, with the amount moved in each growth phase. It hands the result to Cytoscape, to graph formats, to two
matrix formats, and to R as the parameters of a consumer-resource model (CRM), with
[miaSim](https://bioconductor.org/packages/release/bioc/html/miaSim.html) as the example simulator.

It is the sister tool of [grow**net**](https://github.com/crossfeed-bio/crossfeed), which builds
interaction networks from co-cultures, and is built the same way: a thin client with no runtime
dependencies, a local page and a command line, nothing hosted and nothing uploaded.

## Contents

- [Install](#install)
- [Quickstart](#quickstart)
- [What it does](#what-it-does)
- [The outputs](#the-outputs)
- [Consumer-resource models in R](#consumer-resource-models-in-r)
- [The command line](#the-command-line)
- [Guardrails](#guardrails)
- [License](#license)

## Install

With [uv](https://docs.astral.sh/uv/), which fetches a suitable Python by itself:

```bash
uv tool install foodnet
```

or with pipx (`pipx install foodnet`), or from a clone (`pip install -e .`). Each release also carries
`foodnet-<version>-windows.zip`, the whole program in one folder with Python included: unzip it and
double-click `foodnet.exe`. Windows warns about a program few people have run yet; click the small "More
info" link, then "Run anyway".

## Quickstart

```bash
foodnet gui
```

opens the page in your browser. Type taxa in the first box (a species, a strain, a genus or an NCBI taxon
id, one per line), or press Example, and press Get taxon-metabolite network. The result opens with its
buttons (downloads, Send to Cytoscape, the report, the CRM parameters) and an index of its sections; then
come the consumed and produced matrices as an image, the taxa, the arcs and the sources.

![the legend](docs/legend.svg)

## What it does

1. **Taxa to monocultures.** Names resolve to the strains mGrowthDB holds (a species name to all its
   strains). food**net** reads every batch monoculture of those strains that has metabolite measurements.
   A culture-level growth curve counts in a monoculture, since it measures the one strain.
2. **Growth phases.** Exponential growth ends at the first sample where the culture has risen 90% of the
   way from its start to its maximum, or earlier, where the specific growth rate has fallen below a tenth of
   its maximum over two consecutive intervals (a plateau that keeps creeping up, or a late rise that is not
   growth, would otherwise end it late); the stationary phase runs from there to the last metabolite sample.
   A culture without an end of exponential growth (an optical density read without its blank can make one
   look not grown) gives its change over the whole run instead, marked `whole_run`. Below the boxes,
   choose Exponential phase (the default), Stationary phase or Both. A time window in Advanced settings
   replaces the phases, and a second window can be given to metabolites named there (trehalose over the whole
   run, for example), so no compound needs a window of its own. Diauxic shifts are not detected.
3. **A change per phase.** For each replicate and metabolite, the concentration at the end of the phase
   minus the concentration at its start (interpolated between samples), averaged over replicates and then
   over experiments, each experiment counting once. A mean change below the detection limit (0.2 mM, a
   setting; compounds can have their own) is no change, judged on how well the replicates pin the mean
   down: with three or more, a one-sided 90% confidence interval on the mean must lie beyond the limit
   (inside it for no change); with two, both replicates must; one replicate decides nothing. Experiments
   must not contradict each other. Otherwise the value is inconclusive, neither an arc nor a measured
   zero. A metabolite series shorter than 24 h is used, and flagged.
4. **Values from one medium, presence from the others.** With the second box empty, all data are
   considered: the values come from the medium that holds data for the most taxa, and every other medium
   only says whether a compound was produced or consumed (`presence_only` arcs, NA matrix cells). A filled
   second box limits the data to what matches it: study or experiment ids, or a medium name, which gives the
   values from the medium it matches for the most taxa; with ids and no medium, the majority rule runs within
   the ids. "Include supporting evidence outside the
   second box" adds the rest as presence. Because mGrowthDB does not report a medium's composition
   systematically, an experiment whose description says something was added or taken away ("WC plus mucin
   beads", "without glucose", "+Ac"), or whose recorded atmosphere differs, counts as another medium; this
   can be switched off, and experiments can be excluded by id. "Ignore media differences" pools
   every medium; "Report everything as booleans" drops the amounts.
5. **Pooling, with what does not agree reported.** Studies in the value medium are pooled. Experiments
   that disagree on what happened make the value inconclusive, with the caution `conflict`, named in the
   report. The same experiment deposited under two studies is counted once. Because the value medium is
   chosen for the whole search, a taxon's values can change with the other taxa searched: the page names
   the taxa with more data in another medium, and naming a medium in the second box fixes the choice.
6. **Arcs.** A produced arc runs from the taxon to the metabolite, a consumed arc from the metabolite to
   the taxon, and its width is the amount in mM. One arc per study by default; "Merge arcs across studies"
   and "Merge to genus" are advanced settings, as in grow**net**.

**What a taxon produces and consumes alone is not necessarily what it does in a community.** Competition,
cross-feeding, pH and regulation all change it, so the network is a map of capabilities and candidate links,
and a consumer-resource model built from it is a hypothesis to check against the community itself.

**Known limits.** The defaults have not been validated on held-out studies: the checks so far reproduced a
hand-checked matrix with settings chosen for it, and compared directions with a published figure. A net change
hides a compound made and used again within a phase, and one failed sample at a phase end becomes the value.
Growth on an unblanked optical density that starts high can read as no growth; such cultures give their
change over the whole run (whole_run), not a phase. See [docs/METHOD_NOTES.md](docs/METHOD_NOTES.md).

Acid and base forms of one compound are one metabolite (acetic acid and acetate), since mGrowthDB records
both and an HPLC measures one pool. A compound that was never assayed for a taxon is never written as zero.
Every decision and its reason is in [docs/METHOD_NOTES.md](docs/METHOD_NOTES.md).

## The outputs

- **The matrices as an image**, the first section of the result after its buttons and index: consumed and produced side by side (or one
  above the other when they are wide): a number on a gray for a change, white for measured without one,
  pale orange for not assayed, a question mark for inconclusive, a dash for a culture without the phase, an
  open circle for a change seen in another medium, a dot in the corner for a value from one replicate, and a
  dashed frame for a change over the whole run.
  Downloadable as SVG, and in the matrices zip. Hovering over a cell shows its value, its replicates and the
  studies, experiments and medium behind it (also in the downloaded SVG, opened in a browser).
- **Network**: JSON (the canonical format, [schema](schema/metabolite_network.schema.json)) or GraphML.
- **Taxa x metabolites matrix** (CSV): one cell per taxon and metabolite, the mean change in mM, positive
  when produced and negative when consumed. miaSim's efficiency matrix has the opposite signs (positive for
  a resource taken up), so this matrix is not an `E`: the R package builds `E` from the consumed and produced
  matrices, with miaSim's signs.
- **Consumed and produced matrices** (zip): two matrices of non-negative amounts, an evidence matrix for
  each (`measured`, `below_limit`, `whole_run`, `single_replicate`, `seen_elsewhere`, `inconclusive`,
  `no_phase`, `not_grown`, `presence_only`, `not_assayed`; the legend and the help say what each means),
  and a README. In the CRM parameters' stationary phase a second-window compound is `second_window`: its
  change over that window is given once, in the exponential phase. The first header cell of every matrix names the medium the values come
  from.
- **Send to Cytoscape**: the network in a running Cytoscape, in the style of the legend. A downloaded
  GraphML takes the same style from `foodnet style`.
- **Report**: every setting, every arc, and every record left out with its reason.

Every arc also says how long its cultures grew exponentially (`exponential_h`, in hours), in the result
table and in every file.

In every matrix a number is a change beyond the detection limit, 0 is measured without one, and NA is no
value. A change seen only in another medium is NA by default; Advanced settings can write it as TRUE or as
the amount measured there (drawn on a background of its own in the image), and the evidence matrices mark
it `presence_only` either way. With Both, each metabolite has a column per phase.

## Consumer-resource models in R

CRM mode collects growth rates: from the replicates whose metabolites gave the values, else from another
monoculture in the same medium. Get CRM parameters then downloads the matrices, the growth rates, the
initial medium concentrations, each taxon's biomass change over the phase, the hours each value was measured
over (`intervals.csv`) and each amount's lowest and highest replicate (`bounds.csv`), or sends them to R.
With the phase choice Both, the stationary phase (what changed after the end of exponential growth) comes
too, beside the exponential one. Install the companion package once, at the version of your foodnet (the
page and the help print the line for yours, which names the release), and miaSim for the simulations:

```r
install.packages(c("remotes", "BiocManager"))
remotes::install_github("hallucigenia-sparsa/foodnet", subdir = "r", ref = "v0.2.0")
BiocManager::install("miaSim")
```

then:

```r
library(foodnet)
crm <- foodnet_listen()                      # and press Send to R on the page
# NA cells (never assayed, seen only elsewhere, inconclusive) are refused unless you say how: here, as 0
E <- crm_efficiency(crm, na = "zero")                      # miaSim's E, in each taxon's own unit
crm_backcheck(crm, monod_constant = 1, na = "zero")        # each taxon alone against its own monoculture
args <- as_miasim(crm, x0 = crm$biomass_start, monod_constant = 1, missing_resource = 0, na = "zero")
set.seed(1)
tse <- do.call(miaSim::simulateConsumerResource, c(args, list(t_end = 48, t_store = 480)))
abundance <- crm_unscale(crm, SummarizedExperiment::assay(tse))   # back in each growth curve's unit
```

A run is deterministic as `as_miasim()` shapes it. In miaSim 1.18, `simulateConsumerResource` adds a random
immigrant at rate `migration_p` (default 0.01) even with `stochastic = FALSE`, so `as_miasim()` passes
`migration_p = 0`; drift, epochs and external events are off unless `stochastic = TRUE`, and measurement
noise unless `error_variance > 0`. To explore noise, change them in the list before the call, since it already
holds `migration_p` (`args$migration_p <- 0.01; args$stochastic <- TRUE`); to turn it all off again, set
`migration_p = 0, stochastic = FALSE, error_variance = 0`. Leave `norm = FALSE`: relative abundances cannot
be turned back by `crm_unscale()`.

The matrices are measured amounts. miaSim has no uptake rate: a taxon takes up each resource at up to 1 mM
per unit of abundance per hour, so the unit of abundance decides how fast it eats. `crm_efficiency()` and
`as_miasim()` therefore give each taxon a unit of its own (`crm_scale()`), chosen from its biomass change,
growth rate and uptake so that, alone, it grows at its measured rate, gains its measured biomass and makes
its measured by-products in proportion to what it takes up; the same data in another unit give the same
simulation. foodnet measures no Monod constants, and uptake of each resource follows them: choose them, and
check the choice with `crm_backcheck()`, which also compares the time each taxon takes to grow with its
phase. Uptake of each resource is weighted by its share of the measured uptake, so a trace substrate weighs
little, but a taxon grows at its full rate only while its resources are saturating: with real Monod
constants it grows slower than measured, and `crm_backcheck()` shows by how much. `crm_subset()` leaves out a
taxon a simulation cannot use. `as_miasim()` never
lets miaSim draw starting abundances or Monod constants at random, and refuses pooled media and the
stationary phase unless allowed. With the phase choice Both, the CRM is built from the exponential phase,
since a consumer-resource model describes growth; `crm_phase(crm, "stationary")` switches to the other.
See [r/README.md](r/README.md).

## The command line

```bash
foodnet derive --taxa "Escherichia coli LF82" "Bacteroides fragilis" --out network.json --report report.txt
foodnet derive --taxa Roseburia --phase both --format matrices --out roseburia.zip
foodnet derive --taxa Blautia Roseburia --conditions "Wilkins-Chalgren" --crm-mode --crm crm.zip
foodnet validate network.json
```

`foodnet derive --help` lists every option, with examples.

## Guardrails

The same as grow**net**'s: no real data in git (only synthetic fixtures under `tests/fixtures/`), no
runtime dependencies, the network format is a contract checked against its schema, and `make check` runs
lint, the guardrail gate and the tests, in CI and before every commit. See [AGENTS.md](AGENTS.md) and
[CONTRIBUTING.md](CONTRIBUTING.md).

Every arc cites the studies behind it, so attribution resolves at the arc level; see
[docs/DATA_GOVERNANCE.md](docs/DATA_GOVERNANCE.md).

## License

Apache License 2.0; see [LICENSE](LICENSE) and [NOTICE](NOTICE).
