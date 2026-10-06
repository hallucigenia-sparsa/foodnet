"""Export a foodnet network as GraphML, which Cytoscape, igraph, networkx and Gephi read.

The JSON (`foodnet.model.FoodNetwork.to_dict`) is the canonical output; GraphML is the same network for
network tools. Every node keeps its kind (taxon or metabolite) and every arc its direction, phase, evidence,
amount and the studies behind it. It also carries the columns the foodnet Cytoscape style maps (node color,
line style, display width), so a GraphML file given the foodnet style in Cytoscape draws what Send to
Cytoscape draws (grownet learned this the hard way, Karoline, 2026-09-28). A number an arc lacks is left out,
never written as 0.
"""
from __future__ import annotations

from xml.etree import ElementTree as ET

from .model import FoodNetwork

_NS = "http://graphml.graphdrawing.org/xmlns"

# (key id, for, attribute name, attribute type)
_KEYS = [
    ("g_tool", "graph", "tool", "string"),
    ("g_tool_version", "graph", "tool_version", "string"),
    ("g_derived_on", "graph", "derived_on", "string"),
    ("g_derived_at", "graph", "derived_at", "string"),
    ("g_phase", "graph", "phase", "string"),
    ("g_incomplete", "graph", "incomplete", "boolean"),
    ("g_values", "graph", "values", "string"),
    ("n_name", "node", "name", "string"),
    ("n_label", "node", "label", "string"),       # Gephi takes the label from here
    ("n_kind", "node", "kind", "string"),
    ("n_identity", "node", "identity", "string"),
    ("n_taxon_id", "node", "taxon_id", "string"),
    ("n_species", "node", "species", "string"),
    ("n_chebi_id", "node", "chebi_id", "string"),
    ("n_genus", "node", "genus", "string"),
    ("n_color", "node", "color", "string"),
    ("e_direction", "edge", "direction", "string"),
    ("e_phase", "edge", "phase", "string"),
    ("e_evidence", "edge", "evidence", "string"),
    ("e_amount", "edge", "amount", "double"),
    ("e_change", "edge", "change", "double"),
    ("e_sd", "edge", "sd", "double"),
    ("e_n", "edge", "n", "int"),
    ("e_n_experiments", "edge", "n_experiments", "int"),
    ("e_p_value", "edge", "p_value", "double"),
    ("e_q_value", "edge", "q_value", "double"),
    ("e_window_start", "edge", "window_start", "double"),
    ("e_window_end", "edge", "window_end", "double"),
    ("e_exponential_h", "edge", "exponential_h", "double"),
    ("e_medium", "edge", "medium", "string"),
    ("e_study_ids", "edge", "study_ids", "string"),
    ("e_experiments", "edge", "experiments", "string"),
    ("e_cautions", "edge", "cautions", "string"),
    ("e_notes", "edge", "notes", "string"),
    ("e_merged_arcs", "edge", "merged_arcs", "int"),
    ("e_merged_taxa", "edge", "merged_taxa", "string"),
    ("e_line_style", "edge", "line_style", "string"),
    ("e_display_weight", "edge", "display_weight", "double"),
]


def _data(parent, key, value):
    if value is None or value == "":
        return
    d = ET.SubElement(parent, f"{{{_NS}}}data")
    d.set("key", key)
    d.text = str(value)


def to_graphml(net: FoodNetwork, pretty: bool = True) -> str:
    """Serialize the network as a GraphML document (a directed graph)."""
    from .cytoscape import display_weight, genus, line_style, node_colors

    ET.register_namespace("", _NS)
    root = ET.Element(f"{{{_NS}}}graphml")
    for kid, kfor, kname, ktype in _KEYS:
        k = ET.SubElement(root, f"{{{_NS}}}key")
        k.set("id", kid)
        k.set("for", kfor)
        k.set("attr.name", kname)
        k.set("attr.type", ktype)
    graph = ET.SubElement(root, f"{{{_NS}}}graph")
    graph.set("edgedefault", "directed")
    for key in ("tool", "tool_version", "derived_on", "derived_at", "phase", "values"):
        _data(graph, f"g_{key}", net.meta.get(key))
    _data(graph, "g_incomplete", "true" if net.meta.get("incomplete") else "false")
    colors = node_colors(net)
    for node in net.nodes.values():
        n = ET.SubElement(graph, f"{{{_NS}}}node")
        n.set("id", node.id)
        for key, value in (("name", node.name), ("label", node.name or node.id), ("kind", node.kind),
                           ("identity", node.identity), ("taxon_id", node.taxon_id), ("species", node.species),
                           ("chebi_id", node.chebi_id), ("genus", genus(node)), ("color", colors[node.id])):
            _data(n, f"n_{key}", value)
    for i, e in enumerate(net.edges):
        ed = ET.SubElement(graph, f"{{{_NS}}}edge")
        ed.set("id", f"e{i}")
        ed.set("source", e.source)
        ed.set("target", e.target)
        for key, value in (("direction", e.direction), ("phase", e.phase), ("evidence", e.evidence),
                           ("amount", e.amount), ("change", e.change), ("sd", e.sd), ("n", e.n),
                           ("n_experiments", e.n_experiments),
                           ("p_value", e.p_value), ("q_value", e.q_value), ("window_start", e.window_start),
                           ("window_end", e.window_end), ("exponential_h", e.exponential_h), ("medium", e.medium),
                           ("study_ids", " ".join(e.study_ids)), ("experiments", " ".join(e.experiments)),
                           ("cautions", " ".join(e.cautions)), ("notes", "; ".join(e.notes)),
                           ("merged_arcs", e.merged_arcs), ("merged_taxa", " ".join(e.merged_taxa)),
                           ("line_style", line_style(e)), ("display_weight", display_weight(e))):
            _data(ed, f"e_{key}", value)
    if pretty:
        ET.indent(root)
    return '<?xml version="1.0" encoding="UTF-8"?>\n' + ET.tostring(root, encoding="unicode") + "\n"
