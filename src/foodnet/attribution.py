"""Arc-level per-study attribution, as in grownet.

Per-study licenses are respected by citing every study that supports a network at the arc level (each arc
names the studies behind it), rather than bundling (adopted with K. Faust for grownet, 2026-09-14).
"""
from __future__ import annotations

from .model import FoodNetwork


def studies_with_edges(net: FoodNetwork) -> dict:
    """Each study mapped to the arcs it supports."""
    out = {sid: [] for sid in net.studies}
    for e in net.edges:
        for s in e.study_ids:
            out.setdefault(s, []).append(f"{e.source}->{e.target}")
    return out


def render_attribution(net: FoodNetwork) -> str:
    """A human-readable attribution block: studies, licenses, and how many arcs each supports."""
    lines = ["Sources (cited at the arc level; per-study licenses respected):"]
    edges_by_study = studies_with_edges(net)
    for sid, study in sorted(net.studies.items()):
        n = len(edges_by_study.get(sid, []))
        lines.append(f"  {sid}: {study.citation or sid} [{study.license or 'license: see study'}] supports {n} arc(s)")
    return "\n".join(lines)
