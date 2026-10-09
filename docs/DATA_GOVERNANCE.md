# Data governance

## Sources

mGrowthDB is open. foodnet reads batch monocultures (growth curves and metabolite series) from mGrowthDB
through its public API and derives production and consumption itself; mGrowthDB holds no such network.
See https://mgrowthdb.readthedocs.io/en/latest/api.html .

ChEBI (https://www.ebi.ac.uk/chebi, CC BY 4.0) gives each metabolite's formula and charge when growth rates
are on (CRM mode): foodnet asks its public API for the ChEBI ids mGrowthDB records, keeps the answers for the
run only, and names ChEBI and the retrieval date in the CRM parameters (`chemistry_source` in crm.json).

## What is kept on the user's machine

Two things are kept, in the user's cache folder (`~/Library/Caches/foodnet` on macOS,
`%LOCALAPPDATA%\foodnet\Cache` on Windows, `$XDG_CACHE_HOME/foodnet` or `~/.cache/foodnet` elsewhere, or
`FOODNET_CACHE_DIR`), never in the repository (Karoline, 2026-10-09, to query mGrowthDB less):

- **The species list**, for a day: the names, NCBI taxon ids and study ids mGrowthDB holds, which foodnet needs
  to turn typed names into taxa. It is never kept when the crawl could not read every study.
- **What foodnet read of a study** (its experiments, replicates and their series), while the study is unchanged
  and for at most 30 days: mGrowthDB shows only the latest version of a study, and its `uploadedAt` changes with
  every submission (Karoline: "I know that uploadedAt is reliable; it's tied to the submission process"), so a
  search reads the study records from mGrowthDB and uses what is kept only where `uploadedAt` is the same; the
  30 days (Karoline: "30 days") bound what mGrowthDB changes without a submission, such as unit conversions with
  shared metabolite masses. A study without an `uploadedAt` is never kept. The page keeps what it read, study
  records too, in memory for an hour.

Every result says when the species list was read and which studies came from the copy, since when (`data` in
the network JSON and `crm.json`, and the report). `foodnet derive --refresh`, or Read everything from mGrowthDB
again in Advanced settings, reads everything again; deleting the folder does the same. The copies are public
mGrowthDB data under each study's terms (see the attribution below), kept for the user's own searches and not
for redistribution.

## Attribution at the arc level

Per-study licenses are respected by citing every study that supports a network at the arc level: each arc
names the studies behind it, and a network cites all of its supporting studies with their licenses. This
is grownet's rule (adopted with K. Faust, 2026-09-14): it keeps the per-study terms clean and makes clear
exactly which measurements stand behind each arc.

The per-study license is not exposed by the mGrowthDB study endpoint, so arcs carry the study id and url
with the license left empty rather than guessed (the page shows "license: see study").

## Unpublished collaborator data

Any unpublished measurements, media, or genomes shared by a collaborator are used only for the agreed
analysis and never used to train any model. Raw pulled or shared data stays out of version control: see
`.gitignore` and the `no-raw-data` guardrail in `checks/gate.py`. A comparison against an unpublished
figure lives with that figure, not in this repository.

## Outward sign-off

Nothing bearing a collaborator's name or affiliation goes public without their sign-off.
