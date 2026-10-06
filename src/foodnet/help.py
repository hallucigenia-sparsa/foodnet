"""The help page and the About page, as data keyed by the code's own names.

`SETTINGS` has an entry for every key of `foodnet.search.DEFAULTS`, `EDGE_FIELDS` and `NODE_FIELDS` one for
every field of the model, and `tests/test_help.py` requires both, so a new setting or field needs its help
line in the same change (grownet's rule).
"""
from __future__ import annotations

import html

from . import __version__, brand, rbridge
from .brand import COMMAND, NAME, REPOSITORY

ISSUES = f"{REPOSITORY}/issues"
MGROWTHDB = "https://mgrowthdb.gbiomed.kuleuven.be"
GROWNET = "https://github.com/crossfeed-bio/crossfeed"
MIASIM = "https://bioconductor.org/packages/release/bioc/html/miaSim.html"

SETTINGS = {
    "phase": "Which part of growth a change is measured over: the exponential phase (the default), the stationary "
             "phase, or both, each as its own arcs and matrix columns.",
    "window_start": "With a window end, replaces the phases with one window (hours since inoculation).",
    "window_end": "The end of that window. A window that ends after the last metabolite sample uses the last "
                  "sample and marks the value window_beyond_data.",
    "second_window_metabolites": "Metabolites measured over a second time window instead of the phases or the "
                                 "main window, comma separated, by name (acid or base form alike); empty by "
                                 "default. For a compound still being used when the phase ends, such as a "
                                 "slowly consumed sugar.",
    "second_window_start": "The second window's start, in hours since inoculation.",
    "second_window_end": "The second window's end, in hours; empty means until each culture's last sample.",
    "fraction": "Exponential growth ends at the first sample where the culture has risen this share of the way "
                "from its start to its maximum, on the growth curve smoothed by a running median of three (so one "
                "stray point does not move it): 90% by default.",
    "no_growth_factor": "A culture that rose less than this many times did not grow and has no phases (1.5).",
    "detection_limit": "A mean change smaller than this, in either direction, counts as no change (0.2 mM).",
    "judge_spread": "A change must clear the detection limit across the replicates' spread (mean plus and minus "
                    "one standard deviation), not only in its mean; a spread reaching across the limit is "
                    "inconclusive, neither an arc nor a measured zero. On by default; off, the mean alone decides.",
    "compound_limits": "Detection limits of their own for compounds measured at another scale, as name=mM, comma "
                       "separated: thiamine=0.01. Empty by default: every compound takes the detection limit.",
    "ignore_media": "Values from every medium, pooled, instead of only from the value medium.",
    "presence_entries": "How a matrix cell seen only in another medium is written: NA (the cautious default), "
                        "TRUE, or the amount measured there (shown on its own background in the image; not "
                        "comparable with the value medium's amounts). CRM parameters take numbers, so TRUE "
                        "cells reach them as missing.",
    "booleans": "Report 1, 0 or NA instead of amounts; presence in another medium counts as 1.",
    "report_rates": "Collect each taxon's maximum specific growth rate, which a consumer-resource model needs. "
                    "CRM mode switches it on.",
    "rate_method": "How a growth rate is computed: easylinear (as mGrowthDB reports rates) or baranyi.",
    "rate_window": "With easylinear, the number of points in each fitted window (5).",
    "merge_arcs": "One arc per taxon, metabolite, phase and direction across studies, instead of one per study.",
    "min_studies": "Keep only arcs resting on at least this many studies (needs merged arcs above 1).",
    "merge_genera": "One node per genus; a value is the median of its taxa's values, judged against the detection "
                    "limit. Taxa that disagree make it inconclusive.",
    "conditions": "The second box. Study or experiment ids limit the data to them; a medium name chooses the "
                  "value medium (and, alone, limits the data to it). With ids and no medium, the ids' majority "
                  "medium gives the values and their other media presence. Empty: all data, values from the "
                  "medium that holds data for the most taxa, presence from the others.",
    "outside_evidence": "With the second box filled, also take presence-only evidence from every other medium and "
                        "study holding the taxa; the values still come from the second box.",
    "exclude_studies": "Study ids never read.",
    "exclude_experiments": "Experiment ids left out, comma separated; empty by default.",
    "strict_media": "Tell media apart by what an experiment's description or name says was added or taken away "
                    "(\"plus mucin\", \"without glucose\", \"+Ac\") and by a recorded atmosphere, since "
                    "mGrowthDB does not report a medium's composition systematically. On by default; off, only "
                    "the medium names count.",
    "exclude_metabolites": "Metabolites left out by name, comma separated.",
    "include_non_batch": "Command line only: also read chemostat and serial dilution monocultures. Their changes "
                         "are not net changes in the vessel, so this is for exploring, not for values.",
    "spike_factor": "A growth curve with one or two points this many times above both neighbors sets no phase "
                    "boundary (100; 0 switches the check off).",
    "correction": "How the reported q-values are corrected for multiple testing: Benjamini-Hochberg or "
                  "Benjamini-Yekutieli.",
}

NODE_FIELDS = {
    "id": "ncbi:<taxon id> for a strain, the genus and species for a strain without a usable id, genus:<name> "
          "after merging, chebi:<id> for a metabolite (metabolite:<name> without one).",
    "kind": "taxon or metabolite.",
    "name": "The strain's current name in mGrowthDB, or the metabolite's name (acid and base forms are joined "
            "under the base, for example acetic acid under acetate).",
    "identity": "What the id rests on: ncbi, name, genus, chebi or metabolite_name.",
    "taxon_id": "The strain's NCBI taxon id.",
    "species": "Genus and species of a strain's name.",
    "chebi_id": "A metabolite's ChEBI id.",
}

EDGE_FIELDS = {
    "source": "Where the arc starts: the taxon for a produced arc, the metabolite for a consumed one.",
    "target": "Where it ends.",
    "direction": "produced or consumed.",
    "phase": "exponential, stationary, or window.",
    "evidence": "measured (a value from the value medium) or presence_only (seen in another medium, and not in "
                "the value medium).",
    "amount": "mM, the size of the mean net change in this direction. Missing for presence_only arcs and with "
              "booleans; never 0.",
    "change": "mM, the signed mean net change (positive = produced).",
    "sd": "mM, the standard deviation over the experiments' means when several experiments give the value, else "
          "over its replicates.",
    "n": "The replicates behind the value.",
    "n_experiments": "The experiments behind the value: each counts once in the mean, however many replicates it "
                     "has.",
    "p_value": "A one-sample t-test against zero, of the experiments' means (of the replicates when one "
               "experiment gives the value); none when they do not vary. Reported, not used to decide.",
    "q_value": "The p-value corrected for multiple testing over the tests that make the arcs (every study's "
               "values, or the pooled values when arcs are merged).",
    "window_start": "Hours: the mean start of the phase over the replicates.",
    "window_end": "Hours: the mean end of the phase.",
    "exponential_h": "Hours the cultures behind the arc grew exponentially, from their first growth sample to "
                     "the end of exponential growth, averaged over replicates (also with a time window).",
    "medium": "The medium (or media) the arc rests on.",
    "study_ids": "The studies the arc rests on.",
    "experiments": "The mGrowthDB experiments behind it.",
    "cautions": "single_replicate, short_record, window_beyond_data, stationary_not_reached, conflict, "
                "not_detected_in_value_medium, phase_from_other_replicates, coarse_sampling, boundaries_differ, "
                "no_variance (see the legend).",
    "notes": "Remarks in words: what disagreed, what another medium showed.",
    "merged_arcs": "With merged arcs, how many studies the arc joins.",
    "merged_taxa": "With merging to genus, the taxa behind the arc.",
}

SECTIONS = (("what", "What foodnet does"), ("alone", "Alone is not in a community"),
            ("phases", "Growth phases"), ("values", "Values, media and presence"),
            ("matrices", "The two matrix formats"), ("crm", "Consumer-resource models and R"),
            ("cytoscape", "Cytoscape and the downloads"), ("empty", "When a search gives nothing"),
            ("settings", "Every setting"), ("fields", "Every field"), ("cli", "The command line"))


def _e(x) -> str:
    return html.escape(str(x), quote=True)


def _dl(items: dict) -> str:
    return "<dl>" + "".join(f"<dt><code>{_e(k).replace('_', '_<wbr>')}</code></dt><dd>{_e(v)}</dd>"
                            for k, v in items.items()) + "</dl>"


def _default(value) -> str:
    if value is None:
        return "not set"
    if isinstance(value, bool):
        return "on" if value else "off"
    return str(value) if value != "" else "empty"


def render_help(token: str, defaults: dict, example: tuple, job: str = "") -> str:
    toc = "".join(f"<li><a href=\"#{key}\">{_e(title)}</a></li>" for key, title in SECTIONS)
    settings = {k: f"{v} Default: {_default(defaults.get(k))}." for k, v in SETTINGS.items()}
    style_link = f"/foodnet_style.xml?token={_e(token)}"
    return f"""<h2 class="page">Help</h2>
<ul class="toc">{toc}</ul>
<h2 id="what">What foodnet does</h2>
<p>foodnet builds a bipartite network of taxa and metabolites from <a href="{MGROWTHDB}">mGrowthDB</a>: which taxon
produces and which consumes which compound. It reads batch monocultures with metabolite measurements, and an arc's
width is the amount produced or removed. It is the sister tool of <a href="{GROWNET}">grownet</a>, which builds
interaction networks from co-cultures, and works the same way: two boxes, a button, and the results under them.</p>
<p>Try the Example ({_e(", ".join(example))}), or type taxa: a species, a strain, a genus (all its species) or an
NCBI taxon id, one per line.</p>
<h2 id="alone">Alone is not in a community</h2>
<p class="note"><strong>What a taxon produces and consumes on its own is not necessarily what it does in a
community.</strong> Every arc comes from a monoculture: one strain, alone, in one medium, in a batch culture. It
shows what the strain can do under those conditions, not what it will do among other species.</p>
<p>In a community many things change. Partners compete for the same substrates, so a strain may never reach a
compound it uses alone, or may switch to another one it prefers less. Partners also feed it: compounds that are
absent from the medium appear as others produce them, and a strain that only produced a compound alone may take
it up once a partner makes more of it. Its own products, the pH and the densities reached differ, and partners
can switch pathways on or off. A culture that is continuously fed, as in a chemostat, never goes through the
starvation a batch culture ends in, so neither of foodnet's two growth phases need describe it. And a net change
can hide a compound that is produced and consumed at the same time.</p>
<p>So read the network as a map of capabilities and candidate links, not as a prediction. A consumer-resource
model built from these parameters is a hypothesis about the community: check it against measurements of the
community itself before relying on what it predicts.</p>
<h2 id="phases">Growth phases</h2>
<p>What a culture makes or takes up while it grows can differ from what it does once growth has stopped, so every
value belongs to a phase. Exponential growth ends at the first sample at which the culture has risen 90% of the
way from its start to its maximum (an advanced setting), read on the growth curve smoothed by a running median of
three so a single stray point does not move it; the stationary phase runs from there to the last metabolite
sample. A value is the metabolite's concentration at the end of the phase minus its concentration at the start,
interpolated between samples, averaged over replicates, and then over experiments, each experiment counting once.
Diauxic shifts are not detected: a second growth phase counts as stationary. Neither is a slow late rise: a
culture whose plateau keeps creeping up can end its phase late, and when the replicates of one experiment end it
further apart than a sampling interval the value carries the caution boundaries_differ (a time window is then the
safer choice). A boundary placed on fewer than three growth samples carries coarse_sampling. A culture that did
not grow by the no-growth factor has no phases: its cells are no_phase, and a time window gives them values. A
time window in Advanced settings replaces the phases. The growth curve is a per-strain count when the replicate
has one, else the culture's own (cell counts before optical density); a curve that gives no boundary gives way
to the replicate's next one.</p>
<p>A metabolite series shorter than 24 h is still used, and its values carry the caution short_record; the page
warns about it above the result.</p>
<h2 id="values">Values, media and presence</h2>
<p>Values come from one medium: the one that holds data for the most taxa, or the media, experiments or studies
named in the second box. With the second box empty, every other medium only says whether a compound was produced
or consumed: those arcs are presence_only (dashed) and their cells NA in the value matrices. With the second box
filled, study or experiment ids limit the data to them and a medium name chooses the value medium; with ids
and no medium name, the medium holding most taxa within the ids gives the values and their other media
presence. Include supporting evidence outside the second box adds presence from the rest.
Ignore media differences pools every medium.</p>
<p>What counts as one medium: mGrowthDB names a medium but does not report its composition systematically, so an
added or removed compound often appears only in an experiment's description or name ("WC plus mucin beads",
"without glucose", "+Ac"). Such an experiment, or one with another recorded atmosphere, counts as another
medium. A medium name typed in the second box therefore gives the values from the medium it matches for the
most taxa; the others it matches (Wilkins-Chalgren with mucin, for "Wilkins-Chalgren") give presence only, and
the page names them. Tell media apart by their descriptions and atmosphere switches this off; Exclude these
experiments leaves out an experiment by id. The report lists every medium found.
Studies in the value medium are pooled, each experiment counting once; when their experiments disagree on what
happened, the value is inconclusive (NA, no arc), carries the caution conflict, and the report names the
experiments. A change must also clear the detection limit across the replicates' spread: a mean of +0.8 mM with a
standard deviation of 1.0 mM is inconclusive, not production. The same experiment deposited under two studies is
counted once.</p>
<p>Because the value medium is chosen for the whole search, a taxon's values can change with the other taxa
searched with it: when a taxon has more data in another medium, the page says so. Name a medium in the second box
to fix the choice.</p>
<p>Acid and base forms of one compound are one metabolite (acetic acid and acetate), since an HPLC measures the
pool whatever a record calls it. A compound that was not assayed for a taxon is never written as zero.</p>
<h2 id="matrices">The two matrix formats</h2>
<p><strong>Taxa x metabolites (CSV)</strong>: one matrix, rows taxa, columns metabolites, each cell the mean change in
mM, positive when produced and negative when consumed. <strong>Consumed and produced matrices (zip)</strong>: two
matrices of non-negative amounts, with an evidence matrix for each and a README. In both: a number is a change
beyond the detection limit, 0 is measured without one, NA is no value. The evidence matrices say which:
measured, below_limit, seen_elsewhere (0 in the value medium, a change this way in another), inconclusive,
no_phase, presence_only, not_assayed. The first header cell names the medium the values come from. With the
phase choice Both, each metabolite has a column per phase.</p>
<h2 id="crm">Consumer-resource models and R</h2>
<p>CRM mode collects growth rates: from the replicates whose metabolites gave the values, else from another
monoculture in the same medium. Get CRM parameters then downloads the matrices, the rates, the initial medium
concentrations and each taxon's biomass change over the phase with a README, or sends them to R. With Both, the
CRM uses the exponential phase.</p>
<p>The R companion package receives them. Install it once with <code>{_e(rbridge.INSTALL_R)}</code>
({_e(rbridge.INSTALL_TROUBLE)}), then <code>library(foodnet); crm &lt;- foodnet_listen()</code> and press Send to R.
<code>crm_efficiency(crm)</code> builds <a href="{MIASIM}">miaSim</a>'s efficiency matrix from its equations, in a
unit of abundance chosen per taxon (<code>crm_scale</code>) so that, alone, it grows at its measured rate, gains its
measured biomass and makes its measured by-products. Uptake of each resource follows the Monod
constants, which foodnet does not measure: <code>crm_backcheck(crm, monod_constant = 1)</code> simulates each taxon
alone against its monoculture. <code>as_miasim(crm, x0 = crm$biomass_start, monod_constant = 1)</code> gives the
arguments of <code>simulateConsumerResource</code>, and never lets it draw starting abundances or Monod constants at
random.</p>
<h2 id="cytoscape">Cytoscape and the downloads</h2>
<p>Send to Cytoscape puts the network into a running Cytoscape in the style of the legend. A downloaded GraphML can take
the same style: <a href="{style_link}">foodnet_style.xml</a> (File, Import, Styles from File). JSON is the canonical
format; GraphML is the same network for network tools.</p>
<h2 id="empty">When a search gives nothing</h2>
<ul><li>The taxa may have no batch monoculture with metabolites in mGrowthDB: the report lists every record left out
and why.</li>
<li>The second box may match nothing: a medium is matched as text, so a shorter word finds more spellings.</li>
<li>A culture without a growth curve has no phases; a time window does not need one.</li>
<li>Every change may be below the detection limit.</li></ul>
<h2 id="settings">Every setting</h2>{_dl(settings)}
<h2 id="fields">Every field</h2><p>Arcs:</p>{_dl(EDGE_FIELDS)}<p>Nodes:</p>{_dl(NODE_FIELDS)}
<h2 id="cli">The command line</h2>
<pre>{_e(COMMAND)} derive --taxa "Escherichia coli LF82" "Bacteroides fragilis" --out network.json
{_e(COMMAND)} derive --taxa Roseburia --phase both --format matrices --out roseburia.zip
{_e(COMMAND)} derive --taxa Blautia --conditions "Wilkins-Chalgren" --crm-mode --crm crm.zip
{_e(COMMAND)} gui</pre>
<p>Run <code>{_e(COMMAND)} derive --help</code> for every option. Problems and questions:
<a href="{ISSUES}">{_e(ISSUES)}</a>.</p>"""


def render_about() -> str:
    return f"""<h2 class="page">About</h2>
<p>{NAME} {_e(__version__)} builds taxon and metabolite networks from mGrowthDB batch monocultures. It is the sister
tool of grownet and shares its structure: a thin client, no runtime dependencies, nothing hosted, nothing
uploaded. It is developed at the KU Leuven Laboratory of Molecular Bacteriology.</p>
<p>Source and issues: <a href="{REPOSITORY}">{_e(REPOSITORY)}</a>. Data: <a href="{MGROWTHDB}">mGrowthDB</a>, cited
per arc.</p>
<p class="muted">{_e(brand.NAME)} is released under the Apache License 2.0.</p>"""
