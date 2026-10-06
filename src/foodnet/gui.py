"""A local page: type taxa, get which metabolites they produce and consume.

`foodnet gui` starts a small server from the standard library on 127.0.0.1, with a random token in the URL,
and opens the browser. Nothing is hosted, nothing is uploaded, and the page has no JavaScript: the form posts
back and the server returns HTML rendered in Python. This is grownet's page, with the same two input boxes
and the same frame, and foodnet's own settings, results and outputs.

The work is in `foodnet.search.run_query`; rendering is in `render_form` and `render_result`, so both are
testable without a socket.
"""
from __future__ import annotations

import html
import http.server
import json
import secrets
import socketserver
import threading
import time
import urllib.parse
import webbrowser

from . import __version__, brand, matrix, rates, rbridge
from . import help as help_page
from . import selection as selecting
from .attribution import studies_with_edges
from .cytoscape import CytoscapeError, send, style_xml
from .export import to_graphml
from .figure import matrices_svg
from .legend import legend_svg
from .mgrowthdb import MGrowthDBError
from .report import report_text
from .search import DEFAULTS, EXAMPLE, PHASE_LABELS, compound_limits_of, run_query
from .taxonomy import species_index

TITLE = brand.NAME
# one of each kind the first box takes, shown above it
INPUT_EXAMPLES = ("Roseburia intestinalis", "Escherichia coli LF82", "Bacteroides", "536231")
# the network formats the download menu offers (Karoline, 2026-10-04: graphML and JSON "as before", and
# "instead of an adjacency matrix, there will be 2 matrix formats")
FORMATS = (("json", "JSON"), ("graphml", "GraphML"), ("matrix", "Taxa x metabolites matrix (CSV)"),
           ("matrices", "Consumed and produced matrices (zip)"))
EMPTY_HELP = ("foodnet reads batch monocultures with metabolite data. Taxa grown only in co-culture or in a "
              "chemostat, or monocultures without metabolites, give no arcs.")


def _job_suffix(job: str) -> str:
    return f"&job={urllib.parse.quote(job)}" if job else ""


def _esc(x) -> str:
    return html.escape(str(x), quote=True)


def _back(token: str, job: str = "", top: bool = False) -> str:
    href = f"/?token={token}{_job_suffix(job)}" + ("#result" if job else "")
    kind = "bar backtop" if top else "bar"
    return f"<p class=\"{kind}\"><a class=\"btn\" href=\"{_esc(href)}\">Back</a></p>"


HEADING = f"<span class=\"version\">{html.escape(__version__)}</span>"


def _page(body: str, token: str = "", refresh: str = "", job: str = "") -> str:
    """A page in the foodnet style: the header with the mark, the name, the version, Legend, Help, About."""
    t = _esc(token)
    j = _esc(_job_suffix(job))
    icon = urllib.parse.quote(brand.logo_svg(64))
    header = (f"<header><a class=\"brand\" href=\"/?token={t}{j}\"><h1 class=\"brand\">{brand.logo_svg(28)}"
              f"{brand.WORDMARK}</h1></a>{HEADING}"
              f"<nav><a class=\"btn quiet\" href=\"/legend?token={t}{j}\">Legend</a>"
              f"<a class=\"btn quiet\" href=\"/help?token={t}{j}\">Help</a>"
              f"<a class=\"btn quiet\" href=\"/about?token={t}{j}\">About</a></nav></header>")
    reload = f"<meta http-equiv=\"refresh\" content=\"1; url={_esc(refresh)}\">" if refresh else ""
    return ("<!doctype html>\n<html lang=\"en\"><head><meta charset=\"utf-8\">"
            f"<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\"><title>{TITLE}</title>"
            f"{reload}<link rel=\"icon\" href=\"data:image/svg+xml;utf8,{icon}\">"
            f"<style>{brand.CSS}</style></head><body><div class=\"app\">{header}"
            f"<main>{brand.in_prose(body)}</main></div></body></html>\n")


def _checked(flag) -> str:
    return " checked" if flag else ""


def _field(value) -> str:
    return "" if value is None else _esc(value)


def _settings_block(settings: dict) -> str:
    """Every advanced setting, collapsed behind one summary."""
    s = {**DEFAULTS, **settings}
    rate_methods = "".join(f"<option value=\"{m}\"{' selected' if s['rate_method'] == m else ''}>{m}</option>"
                           for m in rates.METHODS)
    presence_options = "".join(
        f"<option value=\"{k}\"{' selected' if s['presence_entries'] == k else ''}>{label}</option>"
        for k, label in (("na", "NA"), ("true", "TRUE"), ("value", "Value from the other medium")))
    corrections = "".join(f"<option value=\"{c}\"{' selected' if s['correction'] == c else ''}>{label}</option>"
                          for c, label in (("bh", "Benjamini-Hochberg"), ("by", "Benjamini-Yekutieli")))
    return f"""<details>
<summary>Advanced settings</summary>
<div class="row"><label>Time window from
  <input name="window_start" type="text" size="5" value="{_field(s['window_start'])}"> h to
  <input name="window_end" type="text" size="5" value="{_field(s['window_end'])}"> h</label>
  <span class="muted">replaces the phase choice above with one window: each metabolite's change from the start
  to the end, in hours since inoculation. Empty by default</span></div>
<div class="row"><label>Another time window for these metabolites
  <input name="second_window_metabolites" type="text" size="24"
   value="{_esc(s['second_window_metabolites'])}"></label>
  <label>from <input name="second_window_start" type="text" size="5" value="{_field(s['second_window_start'])}">
  h to <input name="second_window_end" type="text" size="5" value="{_field(s['second_window_end'])}"> h</label>
  <span class="muted">comma separated metabolite names, measured over this window instead of the phases or the
  window above, for a compound still being used when they end, such as a slowly consumed sugar. An empty end
  means until each culture's last sample. Empty by default</span></div>
<div class="row"><label>Exponential growth ends at
  <input name="fraction" type="text" size="5" value="{_esc(round(s['fraction'] * 100, 6))}"> % of the maximal
  abundance</label>
  <span class="muted">the first sample at which the culture has risen this share of the way from its start to its
  maximum ends the exponential phase and starts the stationary one (90 by default), or earlier, where its growth
  rate has fallen below a tenth of its maximum over two consecutive intervals</span></div>
<div class="row"><label>Detection limit
  <input name="detection_limit" type="text" size="5" value="{_esc(s['detection_limit'])}"> mM</label>
  <span class="muted">a mean change smaller than this, either way, counts as no change; 0.2 by default</span></div>
<div class="row"><label><input type="checkbox" name="judge_confidence" value="1"{_checked(s['judge_confidence'])}>
  Judge changes on a confidence interval</label>
  <span class="muted">three or more replicates: a one-sided 90% confidence interval on the mean beyond the limit;
  two: both beyond it; one decides nothing; and experiments must not contradict each other. Otherwise
  inconclusive (NA, no arc). On by default</span></div>
<div class="row"><label>Evaporation over a run
  <input name="evaporation" type="text" size="5" value="{_esc(round(s['evaporation'] * 100, 6))}"> %</label>
  <span class="muted">a culture whose growth curve shows no growth gives values only for changes beyond this share
  of a compound's level (and beyond the detection limit); 10 by default, more for open plates</span></div>
<div class="row"><label>Detection limits of their own
  <input name="compound_limits" type="text" size="30" value="{_esc(s['compound_limits'])}"></label>
  <span class="muted">for compounds measured at another scale, as name=mM, comma separated: thiamine=0.01</span></div>
<div class="row"><label><input type="checkbox" name="ignore_media" value="1"{_checked(s['ignore_media'])}>
  Ignore media differences</label>
  <span class="muted">values from every medium, pooled. Off by default: values come from the medium that holds data
  for the most taxa (or the media in the second box), and other media only say whether a compound was produced
  or consumed</span></div>
<div class="row"><label><input type="checkbox" name="outside_evidence" value="1"{_checked(s['outside_evidence'])}>
  Include supporting evidence outside the second box</label>
  <span class="muted">off by default: what is entered in the second box limits the search to the data matching it.
  On, every other medium and study holding these taxa adds presence-only evidence (dashed arcs, NA in the
  value matrices); the values still come from the second box. Without anything in the second box, all data
  are considered anyway</span></div>
<div class="row"><label>Entries seen only in another medium
  <select name="presence_entries">{presence_options}</select></label>
  <span class="muted">how a matrix cell is written whose change was seen only in another medium than the values come
  from: NA (default, the cautious choice), TRUE, or the amount measured there, which the image shows on a
  background of its own. The evidence matrices mark these cells presence_only whatever is chosen</span></div>
<div class="row"><label><input type="checkbox" name="booleans" value="1"{_checked(s['booleans'])}>
  Report everything as booleans</label>
  <span class="muted">1 when a taxon produced or consumed a compound (in any medium), 0 when it was measured and
  did not, NA when nothing says either way; no amounts</span></div>
<div class="row"><label><input type="checkbox" name="report_rates" value="1"{_checked(s['report_rates'])}>
  Report growth rates</label>
  <span class="muted">each taxon's maximum specific growth rate in monoculture, as a consumer-resource model takes
  it: from the replicates with metabolite data, else from another monoculture in the same medium. CRM mode
  turns it on</span></div>
<div class="row"><label>Growth rate method
  <select name="rate_method">{rate_methods}</select></label>
  <span class="muted">easylinear (default), the steepest part of the log curve, as mGrowthDB computes the rates it
  reports; or baranyi, a fitted growth model</span></div>
<div class="row"><label>Growth rate window
  <input name="rate_window" type="text" size="5" value="{_esc(s['rate_window'])}"></label>
  <span class="muted">with easylinear: the points in each fitted window (default 5, as mGrowthDB)</span></div>
<div class="row"><label><input type="checkbox" name="merge_arcs" value="1"{_checked(s['merge_arcs'])}>
  Merge arcs across studies</label>
  <span class="muted">one arc per taxon, metabolite, phase and direction, with the value pooled over every study.
  Off by default: one arc per study. The matrices always pool</span></div>
<div class="row"><label><input type="checkbox" name="merge_genera" value="1"{_checked(s['merge_genera'])}>
  Merge to genus</label>
  <span class="muted">one node per genus, each value the median of its taxa's values. Off by default</span></div>
<div class="row"><label>Minimum supporting studies
  <input name="min_studies" type="text" size="5" value="{_esc(s['min_studies'])}"></label>
  <span class="muted">keep arcs resting on at least this many studies; above 1 it needs merged arcs</span></div>
<div class="row"><label>Leave out these metabolites
  <input name="exclude_metabolites" type="text" size="30" value="{_esc(s['exclude_metabolites'])}"></label>
  <span class="muted">comma separated metabolite names; empty by default</span></div>
<div class="row"><label><input type="checkbox" name="strict_media" value="1"{_checked(s['strict_media'])}>
  Tell media apart by their descriptions and atmosphere</label>
  <span class="muted">mGrowthDB names a medium but does not report its composition systematically, so an added or
  removed compound often appears only in the description ("WC plus mucin beads", "without glucose", "+Ac"). On by
  default: such experiments, and ones with another recorded atmosphere, count as another medium. The report lists
  every medium found</span></div>
<div class="row"><label>Exclude these experiments
  <input name="exclude_experiments" type="text" size="30" value="{_esc(s['exclude_experiments'])}"></label>
  <span class="muted">comma separated experiment ids never used; empty by default</span></div>
<div class="row"><label>Exclude these studies
  <input name="exclude_studies" type="text" size="30" value="{_esc(s['exclude_studies'])}"></label>
  <span class="muted">comma separated study ids never read; empty by default</span></div>
<div class="row"><label>Spike limit
  <input name="spike_factor" type="text" size="5" value="{_esc(s['spike_factor'])}"></label>
  <span class="muted">a growth curve with one or two points this many times above both neighbors sets no phase
  boundary; 0 keeps all</span></div>
<div class="row"><label>No-growth factor
  <input name="no_growth_factor" type="text" size="5" value="{_esc(s['no_growth_factor'])}"></label>
  <span class="muted">a culture that rose less than this many times did not grow, so it has no phases; 1.5 by
  default</span></div>
<div class="row"><label>Multiple testing correction
  <select name="correction">{corrections}</select></label>
  <span class="muted">for the q-values reported beside each value (a one-sample t-test of the replicate changes);
  they are reported, never used to decide</span></div>
</details>"""


def _phase_choice(settings: dict) -> str:
    s = {**DEFAULTS, **settings}
    window = s.get("window_start") is not None and s.get("window_end") is not None
    radios = "".join(
        f"<label><input type=\"radio\" name=\"phase\" value=\"{key}\"{_checked(s['phase'] == key)}>{label}</label>"
        for key, label in PHASE_LABELS.items())
    note = (" <span class=\"muted\">(a time window in Advanced settings replaces this)</span>" if window else "")
    return f"<fieldset class=\"phase\"><span class=\"field\">Growth phase</span>{radios}{note}</fieldset>"


def crm_mode_on(settings: dict) -> bool:
    return bool({**DEFAULTS, **(settings or {})}["report_rates"])


def crm_mode(settings: dict, on: bool = True) -> dict:
    """The settings a consumer-resource model needs: a growth rate per taxon (grownet's gLV mode, for CRMs).
    Off puts the setting back to its default."""
    return {**settings, "report_rates": True if on else DEFAULTS["report_rates"]}


CRM_MODE_MESSAGE = ("CRM mode on: growth rates on (Advanced settings). Run the search, then get the CRM "
                    "parameters under the result; press CRM mode again to switch it off.")
CRM_MODE_OFF_MESSAGE = "CRM mode off: growth rates off, which is the default."


def render_form(token: str, entries: str = "", settings: dict | None = None, message: str = "",
                below: str = "", refresh: str = "", job: str = "", conditions: str = "") -> str:
    """The one page: the two boxes, the phase, the settings, and whatever the search produced under them."""
    settings = settings or {}
    note = f"<p class=\"note\">{_esc(message)}</p>" if message else ""
    on = crm_mode_on(settings)
    return _page(f"""{note}<form method="post" action="/run?token={_esc(token)}">
<div class="boxes">
<div class="box">
<label class="field" for="species">Species, strains, genera or NCBI taxon ids</label>
<p class="examples">For example: {" &middot; ".join(_esc(x) for x in INPUT_EXAMPLES)}</p>
<textarea id="species" name="species" rows="5">{_esc(entries)}</textarea>
<p class="hint">One per line (a genus alone stands for all its species), or press Example. Production and
consumption are derived from mGrowthDB batch monocultures on this machine; nothing is uploaded.</p>
</div>
<div class="box">
<label class="field" for="conditions">Media, experiments or studies for the values (optional)</label>
<p class="examples">For example: {" &middot; ".join(_esc(x) for x in selecting.EXAMPLES)}</p>
<textarea id="conditions" name="conditions" rows="5">{_esc(conditions)}</textarea>
<p class="hint">One per line. Empty: all data, values from the medium holding data for the most taxa. Study or
experiment ids limit the data to them; a medium (matched as text, "wilkins" finds every spelling) gives the
values, and without one the ids' majority medium does, their other media giving presence. Advanced settings
can add evidence from outside the box.</p>
</div>
</div>
{_phase_choice(settings)}
<div class="bar"><button class="primary" type="submit">Get taxon-metabolite network</button>
<button type="submit" name="example" value="1">Example</button>
<button type="submit" name="all" value="1">All</button>
<button type="submit" name="crm_mode" value="1" class="switch{' on' if on else ''}" aria-pressed="{str(on).lower()}"
><span class="track"><span class="knob"></span></span>CRM mode</button>
<span class="muted">All reads every batch monoculture with metabolites in mGrowthDB, ignoring the first box.
CRM mode collects growth rates, which a consumer-resource model needs; press it again to switch back.</span></div>
{_settings_block(settings)}
</form>{below}""", token, refresh, job)


def render_token_page(had_token: bool, port: int) -> str:
    what = ("That link carries a token from an earlier run of foodnet, so it no longer opens this one."
            if had_token else
            "This page opens from the link foodnet printed when it started, which carries a one-time token for "
            "this run.")
    return _page(f"""<h2 class="page">foodnet is running on this machine</h2>
<p>{what}</p>
<p class="hint">Where to find the link: the terminal window where foodnet was started prints
<code>foodnet is at http://127.0.0.1:{int(port)}/?token=...</code>. Copy that whole line into the browser.
If the window is gone, stop foodnet there with Ctrl+C and run <code>{_esc(brand.COMMAND)} gui</code> again,
which opens the browser for you.</p>
<p class="hint">The token is what keeps another program on this machine from driving the tool, so this page does
not show it.</p>""")


def render_progress(token: str, job: dict) -> str:
    done, total = job.get("done", 0), job.get("total")
    bar = f"<progress max=\"{total}\" value=\"{done}\"></progress>" if total else "<progress></progress>"
    below = (f"<section class=\"result\" id=\"result\"><h2>Searching</h2>{bar}"
             f"<p class=\"hint\">{_esc(job.get('message', ''))}</p>"
             "<p class=\"hint\">This page updates by itself until the result is ready.</p></section>")
    return render_form(token, "\n".join(job["entries"]), job["settings"], below=below,
                       conditions=(job["settings"] or {}).get("conditions", ""),
                       refresh=f"/?token={token}&job={job['id']}", job=job["id"])


def render_legend(token: str, job: str = "") -> str:
    return _page(_back(token, job, top=True) + f"<div class=\"legend\">{legend_svg()}</div>" + _back(token, job),
                 token, job=job)


def render_help(token: str, job: str = "") -> str:
    return _page(_back(token, job, top=True) + help_page.render_help(token, DEFAULTS, EXAMPLE, job=job)
                 + _back(token, job), token, job=job)


def render_about(token: str, job: str = "") -> str:
    return _page(_back(token, job, top=True) + help_page.render_about() + _back(token, job), token, job=job)


HEADER = ("<tr><th>taxon</th><th>direction</th><th>metabolite</th><th>phase</th><th>mM &plusmn; sd</th>"
          "<th>replicates</th><th>exponential phase (h)</th><th>q</th><th>medium</th><th>remarks</th>"
          "<th>study</th></tr>")


def _amount(e) -> str:
    if e.evidence == "presence_only":
        return ("<span class=\"na\" title=\"seen in another medium: its direction counts, not its size\">"
                "presence only</span>")
    if e.amount is None:
        return "yes"
    return f"{e.amount:.2f}" + ("" if e.sd is None else f" &plusmn; {e.sd:.2f}")


def _hours(value) -> str:
    """The length of the exponential phase, or a quiet dash when no growth curve gave one."""
    if value is None:
        return "<span class=\"na\" title=\"no growth curve gave a boundary\">none</span>"
    return f"{value:.1f}"


def _arc_rows(net, edges) -> str:
    rows = []
    for e in edges:
        pills = [f"caution: {c.replace('_', ' ')}" for c in e.cautions]
        remarks = "".join(f"<span class=\"pill\">{_esc(p)}</span>" for p in pills) + _esc("; ".join(e.notes))
        kind = "made" if e.direction == "produced" else "took"
        q = "" if e.q_value is None else f"{e.q_value:.3g}"
        rows.append(f"<tr><td>{_esc(net.nodes[e.taxon].name)}</td><td class=\"{kind}\">{e.direction}</td>"
                    f"<td>{_esc(net.nodes[e.metabolite].name)}</td><td>{e.phase}</td>"
                    f"<td class=\"num nowrap\">{_amount(e)}</td><td class=\"num\">{_esc(e.n if e.n else '')}</td>"
                    f"<td class=\"num\">{_hours(e.exponential_h)}</td><td class=\"num\">{q}</td>"
                    f"<td>{_esc(e.medium)}</td><td>{remarks}</td>"
                    f"<td>{_esc(' '.join(e.study_ids))}</td></tr>")
    return "".join(rows)


def _web(url: str) -> bool:
    """Whether a study's url (free text from whoever deposited it) is a web address: anything else, such as a
    javascript: url, would run on the page that holds the session token."""
    return urllib.parse.urlparse((url or "").strip()).scheme in ("http", "https")


def _sources(net) -> str:
    if not net.studies:
        return ""
    by_study = studies_with_edges(net)
    items = "".join(
        f"<li>{_esc(sid)}: {_esc(study.citation or sid)} [{_esc(study.license or 'license: see study')}] supports "
        f"{len(by_study.get(sid, []))} arc(s)"
        + (f" &middot; <a href=\"{_esc(study.url)}\">study</a>" if _web(study.url) else "") + "</li>"
        for sid, study in sorted(net.studies.items()))
    return ("<h2>Sources</h2><p class=\"muted\">Cited at the level of each arc; per-study licenses are "
            f"respected.</p><ul class=\"sources\">{items}</ul>")


def _crm_url(token: str, result: dict) -> str:
    port = result.get("port") or ""
    job = result.get("job", "")
    return (f"http://127.0.0.1:{port}/crm.json?token={token}" + (f"&job={job}" if job else "")) if port else ""


def _outputs(token: str, result: dict) -> str:
    t = _esc(token)
    job = _esc(result.get("job", ""))
    tail = f"&amp;job={job}" if job else ""
    options = "".join(f"<option value=\"{k}\">{label}</option>" for k, label in FORMATS)
    download = (f"<form class=\"inline\" method=\"get\" action=\"/download\">"
                f"<input type=\"hidden\" name=\"token\" value=\"{t}\">"
                + (f"<input type=\"hidden\" name=\"job\" value=\"{job}\">" if job else "") +
                "<button class=\"primary\" type=\"submit\">Download network</button> "
                f"<select name=\"format\" aria-label=\"Network format\">{options}</select></form>")
    cytoscape = (f"<form class=\"inline\" method=\"post\" action=\"/cytoscape?token={t}{tail}\">"
                 "<button type=\"submit\">Send to Cytoscape</button></form>")
    report = (f"<details class=\"report\"><summary class=\"btn\">Report</summary>"
              f"<pre>{_esc(report_text(result))}</pre>"
              f"<p><a class=\"btn\" href=\"/report.txt?token={t}{tail}\">Download the report (.txt)</a></p></details>")
    with_rates = bool(result["settings"].get("report_rates"))
    rates_file = (f"<a class=\"btn\" href=\"/rates.csv?token={t}{tail}\">Download the growth rates (.csv)</a>"
                  if result.get("rates") else "")
    crm = (f"<form class=\"inline\" method=\"post\" action=\"/crm?token={t}{tail}\">"
           "<button type=\"submit\">Get CRM parameters</button> "
           "<select name=\"to\" aria-label=\"What to do with the CRM parameters\">"
           "<option value=\"zip\">Download (.zip)</option><option value=\"r\">Send to R</option></select></form>"
           if with_rates and result.get("taxa_nodes") else "")
    tally = matrix.counts(result)
    hint = (f"<p class=\"hint\">The matrices pool every study of the value medium: {tally['measured']} measured "
            f"cell(s), {tally['below_limit']} measured without a change in that direction (0), "
            f"{tally['presence_only']} seen only in another medium and {tally['not_assayed']} not assayed (both "
            "NA). Send to Cytoscape needs Cytoscape running on this machine.</p>")
    r_hint = (f"<p class=\"hint\">Send to R needs an R session waiting for it: install the companion package once "
              f"with <code>{_esc(rbridge.INSTALL_R)}</code>, then run <code>library(foodnet); crm &lt;- "
              f"foodnet_listen()</code> and press this. Without a listener, read the same parameters in R with "
              f"<code>crm &lt;- foodnet_crm(&quot;{_esc(_crm_url(token, result))}&quot;)</code>.</p>"
              if crm else "")
    missing = result.get("without_a_rate") or {}
    missing_note = ""
    if with_rates and missing:
        names = [result.get("names", {}).get(k, k) for k in missing]
        missing_note = (f"<p class=\"hint\">No growth rate for {_esc(', '.join(names[:6]))}"
                        f"{' and more' if len(names) > 6 else ''}: a consumer-resource model needs one from "
                        "elsewhere; the report says why.</p>")
    body = download + cytoscape if result.get("taxa_nodes") else ""
    return (f"<div class=\"bar outputs\">{body}{report}{rates_file}{crm}</div>"
            f"{hint if result.get('taxa_nodes') else ''}{r_hint}{missing_note}")


def _empty_reason(result: dict) -> str:
    if result.get("errors"):
        return "mGrowthDB could not be read completely (see the messages above); try again when it is reachable."
    if not result["resolved"] and not result.get("all"):
        return "None of the entries could be used; each one says why above."
    if not result["studies"]:
        return "mGrowthDB holds these strains, but no study grows them."
    if result["settings"].get("conditions") and not result["value_cultures"]:
        return ("Nothing these taxa were grown in matches the second box (the note above lists what they were "
                "grown in). Change or empty it, or switch on Include supporting evidence outside the second box.")
    if not result["cultures"]:
        return ("The studies holding these taxa have no batch monoculture of them with metabolite data. "
                + EMPTY_HELP)
    if not result["value_cultures"]:
        return ("No monoculture of these taxa matches the second box, so nothing gives values. Empty it, or name "
                "another medium; the report lists the media found.")
    return ("Every change measured stayed below the detection limit, or no culture had a growth curve to place "
            "the phase by; the report says which.")


def _value_note(result: dict) -> str:
    rule = result["value_rule"]
    media = " / ".join(rule["media"]) or "none"
    if rule["rule"] == "all":
        return f"Values from every medium (Ignore media differences): {_esc(', '.join(rule['media']))}."
    if rule["rule"] == "selected":
        if not rule["media"]:
            return "Values: none, since nothing in the second box matched (see the note below)."
        if result["settings"].get("outside_evidence"):
            return f"Values from the second box: {_esc(media)}. Other media and studies give presence only."
        return f"Only data matching the second box: {_esc(media)}."
    n = rule["taxa_per_medium"].get(rule["keys"][0], 0) if rule["keys"] else 0
    if rule["rule"] == "majority_in_scope":
        return (f"Values from {_esc(media)}, the medium holding data for the most taxa ({n}) in the studies or "
                "experiments of the second box. Their other media give presence only.")
    return (f"Values from {_esc(media)}, the medium holding data for the most taxa ({n}). Other media give "
            "presence only.")


def _figure(token: str, result: dict) -> str:
    """The first result (Karoline, 2026-10-04): the consumed and produced matrices as an image, with its
    download."""
    if not result.get("taxa_nodes") or not result.get("metabolite_nodes"):
        return ""
    job = result.get("job", "")
    href = f"/matrices.svg?token={_esc(token)}" + (f"&amp;job={_esc(job)}" if job else "")
    return (f"<h2>Consumed and produced</h2>"
            "<p class=\"muted\">Hover over a cell to see its value, its replicates and the studies, experiments and "
            "medium behind it.</p>"
            f"<div class=\"scroll figure\">{matrices_svg(result)}</div>"
            f"<p class=\"bar\"><a class=\"btn\" href=\"{href}\">Download the image (.svg)</a>"
            "<span class=\"muted\">the same matrices are in the downloads below, as numbers</span></p>")


def _result_section(token: str, result: dict, message: str = "") -> str:
    net = result["network"]
    genera = result.get("genera", {})
    items = []
    for entry, matches in result["resolved"]:
        strains = ", ".join(f"{name} ({tid})" for tid, name in sorted(matches.items()))
        label = f"{_esc(entry)} (genus): {_esc(', '.join(genera[entry]))}<br><span class=\"muted\">" \
                f"{_esc(strains)}</span>" if entry in genera else f"{_esc(entry)}: {_esc(strains)}"
        items.append(f"<li>{label}</li>")
    if result.get("all"):
        items = [f"<li>every batch monoculture with metabolites: {len(result['studies'])} studies read</li>"]
    unresolved = "".join(
        f"<li><strong>{_esc(e)}</strong>: {_esc(result['reasons'].get(e, 'not in mGrowthDB'))}."
        + (f" Did you mean: {_esc(', '.join(result['suggestions'][e]))}?" if result["suggestions"].get(e) else "")
        + "</li>" for e in result["unresolved"])
    unresolved = f"<p class=\"note\">Not used:</p><ul>{unresolved}</ul>" if unresolved else ""
    note = f"<p class=\"note\">{_esc(message)}</p>" if message else ""
    errors = "".join(f"<p class=\"note\">{_esc(e)}</p>" for e in result["errors"])
    warnings = "".join(f"<p class=\"note\">{_esc(w)}</p>" for w in result["warnings"])
    duplicates = "".join(f"<li>{_esc(x)}</li>" for x in result["duplicates"])
    duplicates = (f"<p class=\"muted\">Duplicate deposits, counted once:</p><ul class=\"muted\">{duplicates}</ul>"
                  if duplicates else "")
    outputs = _outputs(token, result)
    edges = sorted(net.edges, key=lambda e: (net.nodes[e.taxon].name, e.phase, e.direction,
                                             net.nodes[e.metabolite].name))
    if edges:
        measured = sum(1 for e in edges if e.evidence == "measured")
        table = (f"<h2>{len(edges)} arc(s): {measured} measured, {len(edges) - measured} presence only</h2>"
                 f"{outputs}<p class=\"muted\">{_value_note(result)}</p>{warnings}{duplicates}"
                 f"<div class=\"scroll\"><table>{HEADER}{_arc_rows(net, edges)}</table></div>")
    else:
        table = (f"<h2>No arcs</h2><p class=\"note\">{_empty_reason(result)} <a href=\"/help?token={_esc(token)}"
                 f"{_esc(_job_suffix(result.get('job', '')))}#empty\">What to try</a>.</p>{warnings}{outputs}")
    studies = ", ".join(result["studies"]) or "none"
    skips = result["skipped"]
    skipped = ""
    if skips:
        shown = "".join(f"<li>{_esc(label)}: {_esc(reason)}</li>" for label, reason in skips[:60]
                        if not reason.startswith("average of the replicates"))
        skipped = (f"<details><summary>{len(skips)} record(s) left out or noted</summary><ul>{shown}</ul>"
                   "<p class=\"muted\">The report lists all of them.</p></details>")
    return (f"<section class=\"result\" id=\"result\">{note}{_figure(token, result)}"
            f"<h2>Taxa</h2><ul>{''.join(items)}</ul>{unresolved}"
            f"<p class=\"muted\">Studies read: {_esc(studies)}</p>{errors}{table}{_sources(net)}{skipped}</section>")


def render_result(token: str, result: dict, message: str = "") -> str:
    settings = result.get("settings") or {}
    return render_form(token, "\n".join(result.get("entries", [])), settings,
                       below=_result_section(token, result, message), job=result.get("job", ""),
                       conditions=settings.get("conditions", ""))


def _number(form, key, cast=float, positive=True):
    raw = (form.get(key, [""])[0] or "").strip()
    if not raw:
        return None
    try:
        value = cast(raw)
    except ValueError:
        return None
    return abs(value) if positive else value


def window_problem(form: dict) -> str:
    """Why a time window typed in Advanced settings cannot be used, or "" when it can (or none was typed)."""
    for prefix, name in (("window", "The time window"), ("second_window", "The second time window")):
        raw_start = (form.get(f"{prefix}_start", [""])[0] or "").strip()
        raw_end = (form.get(f"{prefix}_end", [""])[0] or "").strip()
        if prefix == "window" and not (raw_start or raw_end):
            continue
        if prefix == "second_window" and not raw_end:
            continue                                        # an empty end is "until the last sample"
        try:
            start = float(raw_start or 0)
            end = float(raw_end)
        except ValueError:
            return f"{name} needs numbers in hours; it was {raw_start or '(empty)'} to {raw_end or '(empty)'}."
        if prefix == "window" and not raw_start:
            return f"{name} needs a start as well as an end, in hours."
        if end <= start:
            return f"{name} ends at {end:g} h, which is not after its start at {start:g} h."
    try:
        compound_limits_of({"compound_limits": form.get("compound_limits", [""])[0]})
    except ValueError as e:
        return str(e)[0].upper() + str(e)[1:] + "."
    return ""


def parse_settings(form: dict) -> dict:
    """Settings from the posted form, falling back to the defaults for anything missing or unreadable."""
    s = dict(DEFAULTS)
    phase = form.get("phase", [""])[0]
    if phase in PHASE_LABELS:
        s["phase"] = phase
    start, end = _number(form, "window_start"), _number(form, "window_end")
    if start is not None and end is not None and end > start:
        s["window_start"], s["window_end"] = start, end
    s["second_window_metabolites"] = form.get("second_window_metabolites", [""])[0].strip()
    second_start, second_end = _number(form, "second_window_start"), _number(form, "second_window_end")
    s["second_window_start"] = second_start if second_start is not None else 0.0
    s["second_window_end"] = second_end if second_end is not None and second_end > s["second_window_start"] else None
    evaporation = _number(form, "evaporation")
    if evaporation is not None and 0 <= evaporation < 100:
        s["evaporation"] = evaporation / 100
    fraction = _number(form, "fraction")
    if fraction is not None and 0 < fraction <= 100:
        s["fraction"] = fraction / 100
    for key, cast in (("detection_limit", float), ("spike_factor", float), ("no_growth_factor", float)):
        value = _number(form, key, cast)
        if value is not None:
            s[key] = value
    for key, low in (("rate_window", 2), ("min_studies", 1)):
        value = _number(form, key, int)
        if value is not None:
            s[key] = max(low, value)
    if form.get("rate_method", [""])[0] in rates.METHODS:
        s["rate_method"] = form["rate_method"][0]
    if form.get("presence_entries", [""])[0] in matrix.PRESENCE_ENTRIES:
        s["presence_entries"] = form["presence_entries"][0]
    if form.get("correction", [""])[0] in ("bh", "by"):
        s["correction"] = form["correction"][0]
    for key in ("ignore_media", "booleans", "report_rates", "merge_arcs", "merge_genera", "outside_evidence",
                "strict_media", "judge_confidence"):
        s[key] = bool(form.get(key))
    for key in ("conditions", "exclude_studies", "exclude_experiments", "exclude_metabolites", "compound_limits"):
        s[key] = form.get(key, [""])[0].strip()
    return s


# searches kept for their pages, downloads and reports: the latest ones only
KEPT_JOBS = 20
INDEX_MAX_AGE = 3600
_INDEX_LOCK = threading.Lock()


def prune_jobs(jobs: dict, keep: int = KEPT_JOBS) -> None:
    for old in list(jobs)[:-keep]:
        if jobs[old]["status"] != "running":
            del jobs[old]


class _Handler(http.server.BaseHTTPRequestHandler):
    """Routes: the page, a search and its progress, the downloads, the report, Cytoscape and the CRM
    parameters. Every request carries the session token."""

    token = ""
    client_factory = None
    state: dict = {}
    wait = 1.0

    def log_message(self, *_args):
        pass

    def _send(self, body, content_type: str = "text/html; charset=utf-8", filename: str = "", status: int = 200):
        data = body if isinstance(body, bytes) else body.encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(data)))
        if filename:
            self.send_header("Content-Disposition", f"attachment; filename=\"{filename}\"")
        self.end_headers()
        self.wfile.write(data)

    def _redirect(self, location: str):
        self.send_response(303)
        self.send_header("Location", location)
        self.send_header("Content-Length", "0")
        self.end_headers()

    def _authorized(self, query: dict) -> bool:
        given = query.get("token", [""])[0]
        if secrets.compare_digest(given, self.token):
            return True
        self._send(render_token_page(bool(given), self.server.server_address[1]), status=403)
        return False

    def _job_page(self, job_id: str) -> str:
        job = self.state.setdefault("jobs", {}).get(job_id)
        if job is None:
            return render_form(self.token, message="That search is no longer here; run it again.")
        if job["status"] == "running":
            return render_progress(self.token, job)
        if job["status"] == "failed":
            return render_form(self.token, "\n".join(job["entries"]), job["settings"], job["error"],
                               conditions=(job["settings"] or {}).get("conditions", ""))
        job["result"]["job"] = job["id"]
        job["result"]["port"] = self.server.server_address[1]
        self.state["result"] = job["result"]
        return render_result(self.token, job["result"])

    def _result(self, query: dict):
        """The result a request names by its job, or the last one shown when it names none. A job that is gone
        (pruned, or from before a restart) gives None, never another search's result."""
        wanted = query.get("job", [""])[0]
        if wanted:
            job = self.state.get("jobs", {}).get(wanted)
            result = job["result"] if job and job.get("status") == "done" else None
        else:
            result = self.state.get("result")
        if result is not None:
            result["port"] = self.server.server_address[1]
        return result

    # the routes that serve a search's data: a job they name that is gone is a 404, so an old tab or a pasted
    # crm.json link never receives another search's data
    DATA_ROUTES = ("/download", "/rates.csv", "/crm.zip", "/crm.json", "/matrices.svg", "/report.txt", "/crm",
                   "/cytoscape")

    def _gone(self, path: str, query: dict) -> bool:
        if path in self.DATA_ROUTES and query.get("job", [""])[0] and self._result(query) is None:
            self._send(render_form(self.token, message="That search is no longer here (the page keeps the latest "
                                                       f"{KEPT_JOBS} searches); run it again."), status=404)
            return True
        return False

    def _download(self, fmt: str, query: dict):
        result = self._result(query)
        if not result:
            self._send(render_form(self.token, message="Nothing to download yet."))
        elif fmt == "graphml":
            self._send(to_graphml(result["network"]), "application/xml", f"{TITLE}_network.graphml")
        elif fmt == "matrix":
            self._send(matrix.signed_csv(result["network"], result), "text/csv; charset=utf-8",
                       f"{TITLE}_matrix.csv")
        elif fmt == "matrices":
            self._send(matrix.pair_package(result), "application/zip", f"{TITLE}_matrices.zip")
        else:
            self._send(result["network"].to_json(), "application/json", f"{TITLE}_network.json")

    def do_GET(self):             # noqa: N802 - the name http.server requires
        parsed = urllib.parse.urlparse(self.path)
        query = urllib.parse.parse_qs(parsed.query)
        if not self._authorized(query):
            return
        if self._gone(parsed.path, query):
            return
        if parsed.path == "/":
            job = query.get("job", [""])[0]
            self._send(self._job_page(job) if job else render_form(self.token))
        elif parsed.path in ("/help", "/about", "/legend"):
            job = query.get("job", [""])[0]
            job = job if job in self.state.get("jobs", {}) else ""
            render = {"/help": render_help, "/about": render_about, "/legend": render_legend}[parsed.path]
            self._send(render(self.token, job))
        elif parsed.path == "/download":
            self._download(query.get("format", ["json"])[0], query)
        elif parsed.path in ("/rates.csv", "/crm.zip", "/crm.json"):
            self._crm_files(parsed.path, query)
        elif parsed.path == "/matrices.svg":
            result = self._result(query)
            if not result:
                self._send(render_form(self.token, message="No matrices yet: run a search first."))
            else:
                self._send(matrices_svg(result), "image/svg+xml", f"{TITLE}_matrices.svg")
        elif parsed.path == "/foodnet_style.xml":
            self._send(style_xml(), "application/xml; charset=utf-8", "foodnet_style.xml")
        elif parsed.path == "/report.txt":
            result = self._result(query)
            if not result:
                self._send(render_form(self.token, message="No report yet: run a search first."))
            else:
                self._send(report_text(result), "text/plain; charset=utf-8", f"{TITLE}_report.txt")
        else:
            self.send_error(404, "no such page")

    def _crm_files(self, path: str, query: dict):
        result = self._result(query)
        if not result or not result["settings"].get("report_rates"):
            self._send(render_form(self.token, message="No CRM parameters yet: switch on CRM mode (or Report "
                                                       "growth rates) and run the search again."))
        elif path == "/rates.csv":
            self._send(matrix.rates_csv(result), "text/csv; charset=utf-8", f"{TITLE}_growth_rates.csv")
        elif path == "/crm.json":
            self._send(json.dumps(matrix.crm_payload(result), indent=1), "application/json; charset=utf-8")
        else:
            self._send(matrix.crm_package(result), "application/zip", f"{TITLE}_crm_parameters.zip")

    def _crm(self, query: dict, form: dict) -> None:
        result = self._result(query)
        if not result or not result["settings"].get("report_rates"):
            self._send(render_form(self.token, message="No CRM parameters yet: switch on CRM mode and run the "
                                                       "search again."))
            return
        if form.get("to", ["zip"])[0] != "r":
            self._send(matrix.crm_package(result), "application/zip", f"{TITLE}_crm_parameters.zip")
            return
        payload = matrix.crm_payload(result)
        try:
            answer = rbridge.send(payload)
        except rbridge.RError as e:
            self._send(render_result(self.token, result, message=str(e)))
            return
        note = (f"Sent to R: {answer.get('taxa', 0)} taxa, {answer.get('resources', 0)} resources, "
                f"{answer.get('growth_rates', 0)} growth rate(s). The R session printed what it holds and what to "
                "read before simulating.")
        self._send(render_result(self.token, result, message=note))

    def _to_cytoscape(self, query: dict) -> str:
        result = self._result(query)
        if not result:
            return render_form(self.token, message="Nothing to send yet.")
        try:
            sent = send(result["network"], name=TITLE)
        except CytoscapeError as e:
            return render_result(self.token, result, message=str(e))
        if sent.get("warning"):
            return render_result(self.token, result, message=f"Sent to Cytoscape: network {sent['suid']}; "
                                 f"{sent['warning']}.")
        return render_result(self.token, result, message=f"Sent to Cytoscape: network {sent['suid']}, in the "
                             "foodnet style (as in the legend).")

    def _start(self, entries: list, settings: dict, all_studies: bool = False) -> dict:
        job = {"id": secrets.token_hex(4), "status": "running", "done": 0, "total": None, "message": "Starting",
               "entries": [e.strip() for e in entries if e.strip()], "all": all_studies, "settings": settings,
               "result": None, "error": ""}

        def progress(done, total, message):
            job.update(done=done, total=total, message=message)

        def work():
            try:
                with _INDEX_LOCK:
                    built = self.state.get("index_built", 0)
                    if self.state.get("index") is None or time.monotonic() - built > INDEX_MAX_AGE:
                        fresh = self.client_factory()
                        progress(0, None, "Reading the species list of mGrowthDB (the first search in an hour)")
                        self.state["index"] = species_index(fresh, progress=progress)
                        self.state["client"] = fresh
                        # an index the crawl could not read whole hides taxa: it serves this search (which
                        # reports it) and is built again for the next one
                        self.state["index_built"] = 0 if self.state["index"].failed else time.monotonic()
                    client = self.state["client"]
                job["result"] = run_query(client, entries, settings, self.state["index"], progress=progress,
                                          all_studies=all_studies)
                job["status"] = "done"
            except MGrowthDBError as e:
                job.update(status="failed", error=f"mGrowthDB is not reachable: {e}")
            except Exception as e:             # a bug must reach the page, not only a dead thread
                job.update(status="failed", error=f"The search failed: {type(e).__name__}: {e}")

        jobs = self.state.setdefault("jobs", {})
        jobs[job["id"]] = job
        prune_jobs(jobs)
        job["thread"] = threading.Thread(target=work, daemon=True)
        job["thread"].start()
        return job

    def do_POST(self):            # noqa: N802 - the name http.server requires
        parsed = urllib.parse.urlparse(self.path)
        query = urllib.parse.parse_qs(parsed.query)
        if not self._authorized(query):
            return
        if self._gone(parsed.path, query):
            return
        if parsed.path == "/cytoscape":
            self._send(self._to_cytoscape(query))
            return
        length = int(self.headers.get("Content-Length") or 0)
        form = urllib.parse.parse_qs(self.rfile.read(length).decode("utf-8"))
        if parsed.path == "/crm":
            self._crm(query, form)
            return
        entries = form.get("species", [""])[0].splitlines()
        settings = parse_settings(form)
        if form.get("crm_mode"):
            turning_on = not crm_mode_on(settings)
            settings = crm_mode(settings, turning_on)
            self._send(render_form(self.token, "\n".join(entries), settings,
                                   message=CRM_MODE_MESSAGE if turning_on else CRM_MODE_OFF_MESSAGE,
                                   conditions=settings.get("conditions", "")))
            return
        if form.get("example"):
            self._send(render_form(self.token, "\n".join(EXAMPLE), settings, conditions=settings.get("conditions", "")))
            return
        if form.get("all"):
            job = self._start([], settings, all_studies=True)
            job["thread"].join(self.wait)
            self._redirect(f"/?token={self.token}&job={job['id']}#result")
            return
        if not [e for e in entries if e.strip()]:
            self._send(render_form(self.token, settings=settings, message="Type at least one taxon.",
                                   conditions=settings.get("conditions", "")))
            return
        if window_problem(form):
            # a window that cannot be used is said, not silently dropped (found in the pre-release audit)
            self._send(render_form(self.token, "\n".join(entries), settings, message=window_problem(form),
                                   conditions=settings.get("conditions", "")))
            return
        job = self._start(entries, settings)
        job["thread"].join(self.wait)
        self._redirect(f"/?token={self.token}&job={job['id']}#result")


class _Server(socketserver.ThreadingMixIn, http.server.HTTPServer):
    daemon_threads = True


def serve(port: int = 0, open_browser: bool = True, client_factory=None) -> None:
    """Run the local page until interrupted. port 0 picks a free port."""
    from .mgrowthdb import MGrowthDBClient

    handler = type("FoodnetHandler", (_Handler,), {
        "token": secrets.token_urlsafe(16),
        "client_factory": staticmethod(client_factory or MGrowthDBClient),
        "state": {},
    })
    try:
        server = _Server(("127.0.0.1", port), handler)
    except OSError as e:
        raise SystemExit(f"foodnet gui: cannot use port {port} ({e.strerror}). Pick another with --port, or "
                         "leave it out to use a free one.") from None
    with server as httpd:
        url = f"http://127.0.0.1:{httpd.server_address[1]}/?token={handler.token}"
        print(f"{TITLE} is at {url}\nPress Ctrl+C to stop.", flush=True)
        if open_browser:
            webbrowser.open(url)
        try:
            httpd.serve_forever()
        except KeyboardInterrupt:
            print("\nstopped")
