"""The help page and the About page, as data keyed by the code's own names.

`SETTINGS` has an entry for every key of `foodnet.search.DEFAULTS`, `EDGE_FIELDS` and `NODE_FIELDS` one for
every field of the model, and `tests/test_help.py` requires both, so a new setting or field needs its help
line in the same change (grownet's rule).
"""
from __future__ import annotations

import html

from . import __version__, brand, rbridge
from .brand import COMMAND, NAME, REPOSITORY
from .model import CAUTION_TIERS, CAUTIONS

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
                "from its start to its maximum (90% by default), or earlier where its growth rate has fallen below "
                "a tenth of its maximum over two consecutive intervals.",
    "no_growth_factor": "A culture that rose less than this many times did not grow and has no phases (1.5).",
    "detection_limit": "A mean change smaller than this, in either direction, counts as no change (0.2 mM).",
    "judge_confidence": "With three or more replicates, a change needs a one-sided 90% confidence interval on "
                        "their mean beyond the detection limit (no change: inside it); with two, both replicates "
                        "beyond it on one side (or both inside); one replicate decides nothing. The experiments of "
                        "a value must not contradict each other. Otherwise it is inconclusive, neither an arc nor a "
                        "measured zero. On by default; off, the mean alone decides.",
    "evaporation": "The share of a compound's level that evaporation could change over a run (10%): a culture whose "
                   "growth curve shows no growth counts as metabolically active, and gives values, only with changes "
                   "beyond it. Raise it for open plates or long runs.",
    "compound_limits": "Detection limits of their own for compounds measured at another scale, as name=mM, comma "
                       "separated: thiamine=0.01. Empty by default: every compound takes the detection limit.",
    "ignore_media": "Values from every medium, pooled, instead of only from the value medium.",
    "presence_entries": "How a matrix cell seen only in another medium is written: NA (the cautious default), "
                        "TRUE, or the amount measured there (shown on its own background in the image; not "
                        "comparable with the value medium's amounts). CRM parameters take numbers, so TRUE "
                        "cells reach them as missing.",
    "refresh": "Read everything from mGrowthDB again instead of the copies kept on this machine (the species list, "
               "kept a day; what was read of each study, kept while it is unchanged and at most 30 days); "
               "command line: --refresh.",
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
    "phase": "exponential, stationary, window, or whole_run (a culture without an end of exponential growth: "
             "its change over the whole run).",
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
    "cautions": ", ".join(CAUTIONS) + " (see the legend and Cautions). Ranked: tier 1, the value spans the wrong or an "
                "uneven stretch of time (" + ", ".join(c for c in CAUTIONS if CAUTION_TIERS.get(c) == 1)
                + "); tier 2, whether there is a change, or how large, is less certain ("
                + ", ".join(c for c in CAUTIONS if CAUTION_TIERS.get(c) == 2) + "); tier 3, how the phase was "
                "found, or presence only (" + ", ".join(c for c in CAUTIONS if CAUTION_TIERS.get(c) == 3) + "). The "
                "page names the tier 1 values first, those missing the most first. See Cautions.",
    "notes": "Remarks in words: what disagreed, what another medium showed.",
    "merged_arcs": "With merged arcs, how many studies the arc joins.",
    "merged_taxa": "With merging to genus, the taxa behind the arc.",
}

SECTIONS = (("what", "What foodnet does"), ("alone", "Alone is not in a community"),
            ("phases", "Growth phases"), ("values", "Values, media and presence"), ("cautions", "Cautions"),
            ("matrices", "The two matrix formats"), ("crm", "Consumer-resource models and R"),
            ("miasim", "Simulate with miaSim, step by step"), ("cytoscape", "Cytoscape and the downloads"),
            ("empty", "When a search gives nothing"),
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


def _cautions_of(tier) -> str:
    """The cautions of one tier (None: the ones that explain a missing value), as a list in the legend's words."""
    from .legend import CAUTION_TEXT
    return _dl({c: CAUTION_TEXT[c][:1].upper() + CAUTION_TEXT[c][1:] + "."
                for c in CAUTIONS if CAUTION_TIERS.get(c) == tier})


# the help's account of the cautions (Karoline, 2026-10-09: "add a section that explains the cautions"); every
# caution of foodnet.model.CAUTIONS is listed under its tier, in the legend's words
def _cautions_section() -> str:
    return f"""<h2 id="cautions">Cautions</h2>
<p>A caution is a note on a value, or on why a cell has none. A value with a caution is what its stretch of time
and its replicates show; the caution says what that stretch or those replicates lack. Most values carry
at least one, so cautions are ranked by how much they can change what a value means, from tier 1 (read it first)
to tier 3. A value's cautions are gathered from every replicate and experiment behind it, so one replicate can
bring a caution to a value that three experiments agree on.</p>
<p>Cautions travel with the values: in an arc's <code>cautions</code> field; in <code>cautions.csv</code> in the
downloads, one row for each value or empty cell that carries a caution or a note, the codes separated by spaces
and the notes in words; and in the CRM parameters (in R, <code>crm$caveats$cautions</code>, with the tier of each
code in <code>crm$caveats$caution_tiers</code>). When values carry a tier 1 caution, the page's first warning
names up to eight of them, first those that miss the most of their change after the growth rate fell, then the
other changes by size, and counts the rest. It ranks only cells with a value (produced, consumed or no change):
an empty cell can carry cautions of any tier, which then say why it could not be decided.</p>
<h3>Tier 1: the value spans the wrong or an uneven stretch of time</h3>
<p>Part of the change may fall outside the phase, the value may span the whole run, or its replicates or records
may cover different stretches. The amount can be off, and so can the direction: a compound made and then used
again can show no change, or the other direction, over a stretch that holds both. Read these before using the
value.</p>
<ul>
<li><code>still_changing</code>: the compound kept changing between the two ways of ending growth (see Growth
phases). An exponential value lacks that change, by the mM in its note: too small where the compound kept going
the same way, too large where it turned around. A stationary value starts where the growth rate fell and holds
that change, so an uptake after an earlier rise is understated. With the exponential phase alone (the default),
or where the stationary cell is empty, that change is in no value. A time window over the stretch you need (or
the second time window, for chosen compounds) gives its change in one value.</li>
<li><code>whole_run</code>, <code>boundaries_differ</code>, <code>start_differs</code> and
<code>window_beyond_data</code>: a time window (Advanced settings) inside every replicate's samples gives a value
over one stretch.</li>
<li><code>growth_unclear</code> and <code>growth_unknown</code>: the culture's growth was never confirmed. A time
window, or the second time window for the compounds it names, gives a value too, without this caution, so check
the culture's growth curve in mGrowthDB first.</li>
<li><code>short_record</code>: nothing in the settings helps; the record is short.</li>
</ul>
<p>A time window replaces the phases for the whole search: every culture gets a value over it, also those that
did not grow, and none is marked for that (the warning on cultures that did not grow is not given). The second
time window does the same for the compounds it names, also in a phase search, where that warning still names
such a culture as giving no value. The CRM
back-check cannot judge the Monod constants over it. Compare a window search with the phase search rather than
replace it.</p>
{_cautions_of(1)}
<h3>Tier 2: whether there is a change, or how large, is less certain</h3>
<p>A value with these rests on less: few replicates, identical replicates, experiments that differ in amount, an
experiment left out of the call, a replicate within its series' scatter (how much it jumps from sample to sample
around its trend, checked where the window holds four or more samples) or within what evaporation could do, or a
phase placed from few growth samples or from other replicates. Use it, with less weight. A third replicate helps
most where the value rests on one experiment of two. Identical replicates may be one series deposited twice:
then treat the value as resting on one replicate, which by default decides nothing.</p>{_cautions_of(2)}
<h3>Tier 3: how the phase was found, or presence only</h3>
<p><code>growth_rate_boundary</code> changes nothing by itself: where the compound changed beyond its detection
limit between the two ways of ending growth (on a stationary value, by a quarter of the value's change or more),
the value also carries <code>still_changing</code> (tier 1).
<code>stationary_not_reached</code> leaves a replicate without a stationary value, and the cell empty when none is
left: an uptake that starts after growth was not sampled there, so a missing arc is not evidence of no uptake.
No caution marks a second growth phase (a diauxic shift), which counts as stationary.
<code>not_detected_in_value_medium</code> marks a disagreement between media: the arc's direction comes from
another medium, while the value medium measured no change; weigh that before using it.</p>
{_cautions_of(3)}
<h3>No value</h3>
<p>These say why a cell has no value (NA, and no arc), whatever other cautions it carries; the report and
<code>cautions.csv</code> name the experiments behind it. In the evidence matrices, any other empty cell of an
assayed compound reads <code>no_phase</code>: no stationary phase (<code>stationary_not_reached</code>), or, with
both phases, cultures that reached no end of growth, where <code>cautions.csv</code> says <code>no_phase</code>
too.</p>
{_cautions_of(None)}"""


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
way from its start to its maximum (an advanced setting), or earlier, at the first sample after which the specific
growth rate stays below a tenth of its maximum over two consecutive intervals: a plateau that keeps creeping up,
or a late rise in a qPCR count from lysing cells, would otherwise end the phase days late. The stationary phase
runs from there to the last metabolite sample. A value is the metabolite's concentration at the end of the
phase minus its concentration at the start, interpolated between samples, averaged over replicates, and then
over experiments, each experiment counting once. Diauxic shifts are not detected: a second growth phase counts
as stationary. When the replicates of one experiment end growth further apart than a sampling interval, the
value carries the caution boundaries_differ (a time window is then the safer choice). A boundary placed on fewer
than three growth samples carries coarse_sampling. The growth curve only places the boundary between the
phases: a culture whose curve gives none (one that did not rise by the no-growth factor, as an optical density
read without its blank can make it) gives its change over the whole run instead, marked whole_run, in the
column of the phase asked for, as long as its optical density rose by 0.1 or more, or it used up one compound
and made another beyond the limit (growth_unclear); a culture that did neither gives no value (not_grown), so
drift or a dead inoculum never becomes an arc. A time window in
Advanced settings replaces the phases. The growth curve is a per-strain count when the replicate has one, else
the culture's own (cell counts before optical density); a curve that gives no boundary gives way to the
replicate's next one.</p>
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
experiments. Within an experiment, the replicates must pin the change down: with three or more, a one-sided 90%
confidence interval on their mean must lie beyond the detection limit, so five replicates at -8 mM stay a change
if one sample failed; with two, both must lie beyond it, so -0.2 and +1.9 mM are inconclusive, not production;
one replicate decides nothing. Where the growth rate ended growth more than a sample
before the 90% rule (growth_rate_boundary) and a compound kept changing between the two, its values say
still_changing: the exponential value then lacks part of what the culture took up or made while slowing
(E. coli LF82's glucose). The same experiment deposited under two studies is counted once.</p>
<p>Because the value medium is chosen for the whole search, a taxon's values can change with the other taxa
searched with it: when a taxon has more data in another medium, the page says so. Name a medium in the second box
to fix the choice.</p>
<p>Acid and base forms of one compound are one metabolite (acetic acid and acetate), since an HPLC measures the
pool whatever a record calls it. A compound that was not assayed for a taxon is never written as zero.</p>
{_cautions_section()}
<h2 id="matrices">The two matrix formats</h2>
<p><strong>Taxa x metabolites (CSV)</strong>: one matrix, rows taxa, columns metabolites, each cell the mean change in
mM, positive when produced and negative when consumed (miaSim's efficiency matrix has the opposite signs, so the
R package builds it from the two matrices below, not from this one). <strong>Consumed and produced matrices
(zip)</strong>: two
matrices of non-negative amounts, with an evidence matrix for each and a README. In both: a number is a change
beyond the detection limit, 0 is measured without one, NA is no value. The evidence matrices say which:
measured, below_limit (0), whole_run (a change over the whole run, for a culture without an end of
exponential growth), single_replicate (a value from one replicate), seen_elsewhere (0 in the value medium, a
change this way in another), inconclusive, no_phase (no stationary phase, or no end of growth), not_grown (the
culture did not grow), presence_only, not_assayed. In the CRM parameters' stationary phase, a compound of the
second time window is second_window: its change over that window is given once, in the exponential phase.
The first header cell names the medium the values come from. With the
phase choice Both, each metabolite has a column per phase.</p>
<h2 id="crm">Consumer-resource models and R</h2>
<p>CRM mode collects growth rates: from the replicates whose metabolites gave the values, else from another
monoculture in the same medium. Get CRM parameters then downloads the matrices, the rates, the initial medium
concentrations, each taxon's biomass change over the phase, the hours each value was measured over
(intervals.csv) and bounds on each amount (bounds.csv: its lowest and highest replicate; a 0 from 0 to the
detection limit at least), each resource's formula, charge and degree of reduction from ChEBI (chemistry.csv)
and each taxon's electron balance (electron_balance.csv: the electrons in its measured by-products over those in
what it consumed, with the lowest and highest culture's; a check that changes no value: unmeasured substrates,
unmeasured products and biomass lie outside it, so it has no expected side of 1, and 1 minus it is no biomass
yield; it is withheld when a resource the taxon consumed or produced has no degree of reduction, or when what
it took up is within the assays' detection limits, or when it lies outside its own cultures' range) with a
README, or sends
them to R. For the formulas, foodnet asks ChEBI's public API (www.ebi.ac.uk) in CRM mode.
With Both, the CRM uses the exponential phase and carries the
stationary phase beside it: what changed after the end of exponential growth, where cells may still grow, stop
or die (its biomass change says which); in R,
<code>crm_phase(crm, "stationary")</code> switches to it.</p>
<p>The R companion package receives them. Install it once with <code>{_e(rbridge.INSTALL_R)}</code>
({_e(rbridge.INSTALL_TROUBLE)}), then <code>library(foodnet); crm &lt;- foodnet_listen()</code> and press Send to R.
The steps of a simulation are <a href="#miasim">below</a>. <code>crm_efficiency(crm)</code> builds
<a href="{MIASIM}">miaSim</a>'s efficiency matrix from its equations, in a
unit of abundance chosen per taxon (<code>crm_scale</code>) so that, alone, it grows at its measured rate, gains its
measured biomass and makes its measured by-products. Uptake of each resource follows the Monod
constants, which foodnet does not measure: <code>crm_backcheck(crm, monod_constant = 1)</code> simulates each taxon
alone against its monoculture. <code>as_miasim(crm, x0 = crm$biomass_start, monod_constant = 1)</code> gives the
arguments of <code>simulateConsumerResource</code>, and never lets it draw starting abundances or Monod constants at
random.</p>
<h2 id="miasim">Simulate with miaSim, step by step</h2>
<p>Press <strong>CRM example</strong> to follow these steps with five gut species in Wilkins-Chalgren (studies 2, 4,
7 and 9, 0 to 48 h, trehalose over the whole run); with your own taxa, switch CRM mode on instead. The example
uses a time window, so every taxon has a growth rate, a biomass gain and an uptake and the steps run as they are;
with the exponential phase (the default) the same taxa run too, and only there does the back-check (step 4) judge
the Monod constants.</p>
<ol>
<li>Install the R package and <a href="{MIASIM}">miaSim</a> once: <code>install.packages(c("remotes",
"BiocManager"))</code>, <code>{_e(rbridge.INSTALL_R)}</code>, <code>BiocManager::install("miaSim")</code>
({_e(rbridge.INSTALL_TROUBLE)}).</li>
<li>In R: <code>library(foodnet); crm &lt;- foodnet_listen()</code>. On this page press Get taxon-metabolite network,
then Get CRM parameters with Send to R. (Without a listener: <code>crm &lt;- foodnet_crm(url)</code>, with the
address shown under the button.)</li>
<li><code>print(crm)</code>, and read what it says first: the values whose cautions change what they mean, and any
taxon without a growth rate, a biomass gain or an uptake. A simulation cannot use such a taxon; leave it out with
<code>crm &lt;- crm_subset(crm, taxa = setdiff(crm$taxa, c("its name")))</code> (the error names it, if you go
on).</li>
<li>Choose Monod constants, which foodnet does not measure, and check them: <code>crm_backcheck(crm, monod_constant
= 1, na = "zero")</code> simulates each taxon alone and compares it with its monoculture. Read two rows per taxon:
"biomass" (a ratio near 1: it gains what it gained over the phase) and "hours to grow" (a ratio near 1: it takes
as long as it did; below 1, faster; NA: it never gains all of it, too slowly or levelling off just short). These
judge the constants over the exponential phase, whose length is the time each taxon grew. Over a time window, as
in the CRM example, they do not: the measured hours are the window's length, and the biomass ratio is near 1
while the constants lie well below the resource concentrations, so <code>crm_backcheck()</code> warns. To choose
them, search the same taxa with the exponential phase (switch the window off in Advanced settings) and check
there. A ratio that stays put as the constants fall toward 0 is not theirs to fix: the taxon then grows as fast
as the data let it, and what is left lies in the data. Its growth rate may be below its phase's own mean rate
(<code>print(crm)</code> says so; <code>growth = "phase_floor"</code>, given to <code>crm_backcheck()</code> and
<code>as_miasim()</code> alike, uses the phase's rate there), or part of its uptake may come after the growth rate
fell (still_changing, which the warnings name). It takes <code>crm</code> itself and builds the efficiency matrix;
<code>na = "zero"</code> counts NA cells as 0 and says how many of each kind.</li>
<li>Simulate the community with the parameters and the arguments you checked in step 4 (the same
<code>monod_constant</code>, and <code>growth</code> if you set it): <code>args &lt;- as_miasim(crm, x0 =
crm$biomass_start, monod_constant = 1, missing_resource = 0, na = "zero")</code>, then <code>tse &lt;-
do.call(miaSim::simulateConsumerResource, c(args, list(t_end = 48, t_store = 480)))</code>.</li>
<li>Back in each growth curve's unit: <code>abundance &lt;- crm_unscale(crm, SummarizedExperiment::assay(tse),
args = args)</code>, which takes the unit the simulation used from <code>args</code>, plotted over hours with
<code>matplot(SummarizedExperiment::colData(tse)$time, t(log10(abundance)), type = "l")</code>.</li>
</ol>
<p>The run is deterministic as <code>as_miasim()</code> shapes it: miaSim adds random immigrants at
<code>migration_p</code> even without <code>stochastic</code>, so it passes <code>migration_p = 0</code>. To explore
noise, change <code>args</code> before the call (<code>args$migration_p &lt;- 0.01; args$stochastic &lt;-
TRUE</code>); to turn it off again, <code>args$migration_p &lt;- 0; args$stochastic &lt;- FALSE;
args$error_variance &lt;- 0</code>.</p>
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
{_e(COMMAND)} derive --taxa Escherichia coli LF82, Bacteroides fragilis --out network.json
{_e(COMMAND)} derive --taxa Roseburia --phase both --format matrices --out roseburia.zip
{_e(COMMAND)} derive --taxa Blautia --conditions "Wilkins-Chalgren" --crm-mode --crm crm.zip
{_e(COMMAND)} gui</pre>
<p>The first search of a day reads mGrowthDB's species list from every study, which takes a few minutes when
mGrowthDB is slow (the progress text then says how many requests were slow or retried); the list is kept for a
day in your cache folder, and what foodnet read of a study is kept there too, used again while mGrowthDB shows
the study unchanged (the same <code>uploadedAt</code>) and at most 30 days. Read everything from mGrowthDB again
(Advanced settings), or <code>{_e(COMMAND)} derive --refresh</code>, uses neither. The folder:
<code>~/Library/Caches/foodnet</code> on macOS, <code>%LOCALAPPDATA%\\foodnet\\Cache</code> on Windows,
<code>$XDG_CACHE_HOME/foodnet</code> or <code>~/.cache/foodnet</code> elsewhere, or <code>FOODNET_CACHE_DIR</code>;
deleting it reads everything
again.</p>
<p>Run <code>{_e(COMMAND)} derive --help</code> for every option. Problems and questions:
<a href="{ISSUES}">{_e(ISSUES)}</a>.</p>"""


# the released versions, newest first (Karoline, 2026-10-09: "include the release history in the About page");
# tests/test_help.py checks them against CHANGELOG.md
RELEASES = (
    ("0.4.0", "2026-10-09", "faster mGrowthDB access, kept species list and studies, --refresh"),
    ("0.3.0", "2026-10-08", "chemistry and electron balance in CRM mode"),
    ("0.2.0", "2026-10-07", "intervals, bounds, stationary phase and ranked cautions in the CRM parameters"),
    ("0.1.0", "2026-10-06", "first release"),
)


def _changelog_ref() -> str:
    """The changelog of this version: its tag for a release, main for a development version (an audit: a link to
    main showed changes the installed version does not have)."""
    return "main" if "dev" in __version__ else f"v{__version__}"


def render_about() -> str:
    return f"""<h2 class="page">About</h2>
<p>{NAME} {_e(__version__)} builds taxon and metabolite networks from mGrowthDB batch monocultures. It is the sister
tool of grownet and shares its structure: a thin client, no runtime dependencies, nothing hosted, nothing
uploaded. It is developed at the KU Leuven Laboratory of Molecular Bacteriology.</p>
<p>Source and issues: <a href="{REPOSITORY}">{_e(REPOSITORY)}</a>. Data: <a href="{MGROWTHDB}">mGrowthDB</a>, cited
per arc; in CRM mode, the metabolites' formulas and charges from <a href="https://www.ebi.ac.uk/chebi">ChEBI</a>
(CC BY 4.0), which foodnet then contacts as well.</p>
<h3>Release history</h3>
<ul>{"".join(f"<li><b>{_e(v)}</b> ({_e(d)}): {_e(what)}</li>" for v, d, what in RELEASES)}</ul>
<p>Every change, version by version: <a href="{REPOSITORY}/blob/{_changelog_ref()}/CHANGELOG.md">CHANGELOG.md</a>.</p>
<p class="muted">{_e(brand.NAME)} is released under the Apache License 2.0.</p>"""
