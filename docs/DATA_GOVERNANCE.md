# Data governance

## Sources

mGrowthDB is open. foodnet reads batch monocultures (growth curves and metabolite series) from mGrowthDB
through its public API and derives production and consumption itself; mGrowthDB holds no such network.
See https://mgrowthdb.readthedocs.io/en/latest/api.html .

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
