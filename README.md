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
id, one per line), or press Example, and press Get taxon-metabolite network. The consumed and produced
matrices appear first under the settings, then the taxa and the arcs with the downloads.

![the legend](docs/legend.svg)

## What it does

1. **Taxa to monocultures.** Names resolve to the strains mGrowthDB holds (a species name to all its
   strains). food**net** reads every batch monoculture of those strains that has metabolite measurements.
   A culture-level growth curve counts in a monoculture, since it measures the one strain.
2. **Growth phases.** Exponential growth ends at the first sample where the culture has risen 90% of the
   way from its start to its maximum, read on the growth curve smoothed by a running median of three; the
   stationary phase runs from there to the last metabolite sample. Below the boxes,
   choose Exponential phase (the default), Stationary phase or Both. A time window in Advanced settings
   replaces the phases, and a second window can be given to metabolites named there (trehalose over the whole
   run, for example), so no compound needs a window of its own. Diauxic shifts are not detected.
3. **A change per phase.** For each replicate and metabolite, the concentration at the end of the phase
   minus the concentration at its start (interpolated between samples), averaged over replicates and then
   over experiments, each experiment counting once. A mean change below the detection limit (0.2 mM, a
   setting; compounds can have their own) is no change, and a change must clear the limit across the
   replicates' spread: a spread reaching across it is inconclusive, neither an arc nor a measured zero. A
   metabolite series shorter than 24 h is used, and flagged.
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

Acid and base forms of one compound are one metabolite (acetic acid and acetate), since mGrowthDB records
both and an HPLC measures one pool. A compound that was never assayed for a taxon is never written as zero.
Every decision and its reason is in [docs/METHOD_NOTES.md](docs/METHOD_NOTES.md).

## The outputs

- **The matrices as an image**, shown first under the settings: consumed and produced side by side (or one
  above the other when they are wide): a number on a gray for a change, white for measured without one,
  pale orange for not assayed, a question mark for inconclusive, a dash for a culture without the phase, an
  open circle for a change seen in another medium.
  Downloadable as SVG, and in the matrices zip. Hovering over a cell shows its value, its replicates and the
  studies, experiments and medium behind it (also in the downloaded SVG, opened in a browser).
- **Network**: JSON (the canonical format, [schema](schema/metabolite_network.schema.json)) or GraphML.
- **Taxa x metabolites matrix** (CSV): one cell per taxon and metabolite, the mean change in mM, positive
  when produced and negative when consumed.
- **Consumed and produced matrices** (zip): two matrices of non-negative amounts, an evidence matrix for
  each (`measured`, `below_limit`, `seen_elsewhere`, `inconclusive`, `no_phase`, `presence_only`,
  `not_assayed`), and a README. The first header cell of every matrix names the medium the values come
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
initial medium concentrations and each taxon's biomass change over the phase, or sends them to R. Install
the companion package once:

```r
remotes::install_github("hallucigenia-sparsa/foodnet", subdir = "r")
```

then:

```r
library(foodnet)
crm <- foodnet_listen()                      # and press Send to R on the page
E <- crm_efficiency(crm)                     # miaSim's yields and by-products per unit of growth
crm_backcheck(crm, monod_constant = 1)       # each taxon alone against its own monoculture
args <- as_miasim(crm, x0 = crm$biomass_start, monod_constant = 1)
set.seed(1)
tse <- do.call(miaSim::simulateConsumerResource, c(args, list(t_end = 48, t_store = 480)))
```

The matrices are measured amounts. miaSim reads a positive entry of its efficiency matrix as a yield (the
biomass made per mM taken up; uptake itself follows the Monod constants) and a negative one as a by-product
per unit of growth, so `crm_efficiency()` builds E from the amounts, the biomass changes and the growth
rates such that a taxon alone gains its measured biomass and makes its measured by-products. foodnet
measures no Monod constants: choose them, and check the choice with `crm_backcheck()`. `as_miasim()` never
lets miaSim draw starting abundances or Monod constants at random, and refuses pooled media and the
stationary phase unless allowed. With the phase choice Both, the CRM parameters use the exponential phase,
since a consumer-resource model describes growth. See [r/README.md](r/README.md).

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
