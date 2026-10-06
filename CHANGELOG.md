# Changelog

All notable changes to foodnet. The format follows Keep a Changelog, and the version numbers follow
semantic versioning.

## [0.1.0] (unreleased)

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
- A filled second box limits the search: study or experiment ids set the scope, a medium name the value
  medium, and with ids alone the majority medium within them gives the values and their other media
  presence. The advanced option "Include supporting evidence outside the second box" adds the rest as
  presence. An empty box considers all data.
- A second time window for metabolites named in Advanced settings (an empty end means the last sample), so
  a compound such as trehalose can be measured over the whole run while the rest use the phases or a window.
- Strains are named by their current name in mGrowthDB, as in grownet.
- A second box that matches none of the taxa's monocultures says so and lists the media they were grown in.
- No metabolite is left out by default; "Leave out these metabolites" is an advanced setting.
