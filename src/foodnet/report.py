"""The report of a search, as plain text: every setting, every arc, and every reason something gave nothing.

The page shows it under Report and serves it as report.txt; the command line writes it with --report. It
holds what a reader needs to judge the network without the page: which data gave values and which only
presence, the duplicates counted once, the conflicts pooled, and why each left-out record was left out.
"""
from __future__ import annotations

from collections import Counter

from .attribution import render_attribution
from .matrix import conflicts, counts, presence_lines


def _value(value) -> str:
    if value is None:
        return "not set"
    if isinstance(value, bool):
        return "yes" if value else "no"
    return str(value)


def _arc_line(net, e) -> str:
    taxon = net.nodes[e.taxon].name
    met = net.nodes[e.metabolite].name
    if e.amount is not None:
        size = f"{e.amount:.2f} mM" + (f" (sd {e.sd:.2f}, n {e.n})" if e.sd is not None else f" (n {e.n})")
    elif e.evidence == "presence_only":
        size = "presence only"
    else:
        size = f"yes (n {e.n})"
    extra = (f"; cautions: {', '.join(e.cautions)}" if e.cautions else "") + \
        (f"; {'; '.join(e.notes)}" if e.notes else "")
    grew = "" if e.exponential_h is None else f"; exponential phase {e.exponential_h:.1f} h"
    return (f"  {taxon} {e.direction} {met}, {e.phase} phase: {size}{grew}; {e.medium} "
            f"[{' '.join(e.study_ids)}]{extra}")


def report_text(result: dict) -> str:
    net = result["network"]
    meta = net.meta
    s = result["settings"]
    rule = result["value_rule"]
    lines = [f"foodnet {meta.get('tool_version', '')} report, derived {meta.get('derived_at', '')}",
             f"Source: {meta.get('source_db', 'mGrowthDB')}", ""]
    if result.get("all"):
        lines.append(f"Query: every batch monoculture with metabolites in mGrowthDB ({len(result['studies'])} studies)")
    else:
        lines.append("Query: " + ", ".join(result["entries"]))
        for entry, matches in result["resolved"]:
            lines.append(f"  {entry}: " + ", ".join(f"{n} ({t})" for t, n in sorted(matches.items())))
        for entry in result["unresolved"]:
            lines.append(f"  {entry}: not used, {result['reasons'].get(entry, 'not in mGrowthDB')}")
    lines += ["", "Settings:"] + [f"  {k}: {_value(v)}" for k, v in sorted(s.items())]
    lines += ["", f"Studies read: {', '.join(result['studies']) or 'none'}",
              f"Replicates with metabolite data: {result['cultures']}, of which {result['value_cultures']} give values"]
    if rule["rule"] in ("majority", "majority_in_scope"):
        where = " in the studies or experiments of the second box" if rule["rule"] == "majority_in_scope" else ""
        lines.append(f"Values from the medium holding data for the most taxa{where}: "
                     + " / ".join(rule["media"] or ["none"]))
    elif rule["rule"] == "selected":
        lines.append("Values from the second box: " + ", ".join(rule["media"] or ["nothing matched"]))
    else:
        lines.append("Values from every medium (Ignore media differences)")
    for key, n in sorted(rule["taxa_per_medium"].items(), key=lambda kv: -kv[1]):
        lines.append(f"  {key}: {n} taxon(s)")
    if result["warnings"]:
        lines += ["", "Warnings:"] + [f"  * {w}" for w in result["warnings"]]
    if result["errors"]:
        lines += ["", "Errors:"] + [f"  * {e}" for e in result["errors"]]
    if result["duplicates"]:
        lines += ["", "Duplicate deposits, counted once:"] + [f"  * {x}" for x in result["duplicates"]]
    found = conflicts(result)
    if found:
        lines += ["", "Values that pool experiments that disagree:"] + [f"  * {x}" for x in found]
    tally = counts(result)
    lines += ["", f"Matrix cells (consumed and produced together): {tally['measured']} measured, "
              f"{tally['below_limit']} measured below the limit or the other way, {tally['presence_only']} "
              f"presence only, {tally['not_assayed']} not assayed"]
    measured = [e for e in net.edges if e.evidence == "measured"]
    presence = [e for e in net.edges if e.evidence == "presence_only"]
    lines += ["", f"{len(net.edges)} arc(s): {len(measured)} measured, {len(presence)} presence only"]
    lines += [_arc_line(net, e) for e in sorted(measured + presence,
                                                  key=lambda e: (net.nodes[e.taxon].name, e.phase, e.direction,
                                                                 net.nodes[e.metabolite].name))]
    other = presence_lines(result)
    if other:
        lines += ["", "Seen in other media (presence only):"] + [f"  * {x}" for x in other]
    if s.get("report_rates"):
        lines += ["", "Growth rates (1/h):"]
        for tid, r in sorted(result["rates"].items(), key=lambda kv: net.nodes[kv[0]].name):
            lines.append(f"  {net.nodes[tid].name}: {r['rate']:.3f} (median of {r['n']}, from {r['source']}, "
                         f"{r['method']})")
        for tid, why in result["without_a_rate"].items():
            lines.append(f"  {net.nodes[tid].name if tid in net.nodes else tid}: none, {why}")
    if result["initial"]:
        lines += ["", "Initial concentrations in the value medium (mM, mean of the first samples):"]
        for _, v in sorted(result["initial"].items(), key=lambda kv: kv[1]["name"]):
            lines.append(f"  {v['name']}: {v['mean']:.3g} (range {v['min']:.3g} to {v['max']:.3g}, n {v['n']})")
    if result["skipped"]:
        lines += ["", f"Left out or noted ({len(result['skipped'])}):"]
        reasons = Counter(r for _, r in result["skipped"])
        for label, reason in result["skipped"]:
            if reasons[reason] > 3 and reason.startswith("average of the replicates"):
                continue
            lines.append(f"  * {label}: {reason}")
        averages = sum(n for r, n in reasons.items() if r.startswith("average of the replicates"))
        if averages > 3:
            lines.append(f"  * {averages} average replicate(s) left out (average of the replicates, not an "
                         "independent replicate)")
    lines += ["", render_attribution(net)]
    return "\n".join(lines) + "\n"
