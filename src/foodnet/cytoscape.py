"""Send a network into a running Cytoscape, with the style the legend describes.

Cytoscape ships CyREST, a REST API on `http://localhost:1234/v1`, so a plain POST puts a derived network
into the open session: no file, no dependency, nothing hosted. This is grownet's route, with the CyREST
lessons it learned kept (a style is applied with GET; an existing style is updated in place; a number an
arc lacks is left out, not sent as null, and its column declared).

What the style draws:

  * taxa as circles colored by genus, metabolites as rounded squares in a pale neutral, so the two kinds
    of node read apart at a glance and the picture is visibly bipartite;
  * produced arcs blue and consumed arcs amber, from the taxon to the metabolite and from the metabolite to
    the taxon, all ending in the same arrowhead, so the color and the direction say the same thing twice;
  * width by `display_weight`, the amount in mM, so a thicker arc moved more of the compound; an arc with
    no amount (presence only, or booleans) is drawn at a fixed width;
  * long dashes for presence_only arcs (seen in another medium, so its size is unknown) and dots for a
    single replicate, combined when both;
  * stationary-phase arcs fainter than exponential ones, so a "Both" search shows the two phases in one
    picture without hiding either;
  * no edge labels; every number is a column, and Label can be mapped to `amount` in the Style tab.

Everything is sent to the configured localhost port and nowhere else.
"""
from __future__ import annotations

import http.client
import json
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter

from . import brand
from .model import FoodNetwork, genus_name

PORT = 1234
STYLE_NAME = brand.NAME
FIXED_WIDTH = 2.5             # an arc without an amount
WIDE_AT_MM = 20.0             # the amount drawn at the widest


def base_url(port: int = PORT) -> str:
    return f"http://127.0.0.1:{port}/v1"


def line_style(edge) -> str:
    """The dash pattern of an arc: one column, because Cytoscape maps one column per property."""
    presence = edge.evidence == "presence_only"
    single = "single_replicate" in edge.cautions
    if presence and single:
        return "DASH_DOT"
    if presence:
        return "LONG_DASH"
    if single:
        return "DOT"
    return "SOLID"


def display_weight(edge) -> float:
    """The width column: the amount in mM, or a fixed width for an arc that has none."""
    return float(edge.amount) if edge.amount is not None else FIXED_WIDTH


def genus(node) -> str:
    return genus_name(node.name or node.species) if node.kind == "taxon" else ""


def node_colors(net: FoodNetwork) -> dict:
    """node id -> fill color: each genus its own (the most common first), every metabolite the same pale."""
    counts = Counter(genus(n) for n in net.taxa())
    ranked = sorted(counts, key=lambda g: (-counts[g], g))
    by_genus = {g: brand.GENUS_COLORS[i % len(brand.GENUS_COLORS)] for i, g in enumerate(ranked)}
    return {n.id: by_genus[genus(n)] if n.kind == "taxon" else brand.METABOLITE_NODE for n in net.nodes.values()}


def network_json(net: FoodNetwork, name: str = "foodnet") -> dict:
    """The network as Cytoscape.js JSON: every node and arc attribute becomes a column."""
    colors = node_colors(net)
    nodes = [{"data": {k: v for k, v in {
        "id": n.id, "name": n.name or n.id, "kind": n.kind, "identity": n.identity, "taxon_id": n.taxon_id,
        "species": n.species, "chebi_id": n.chebi_id, "genus": genus(n), "color": colors[n.id]}.items()}}
        for n in net.nodes.values()]
    edges = []
    for e in net.edges:
        data = {"source": e.source, "target": e.target,
                # the interaction column decides whether Cytoscape merges parallel arcs, so it carries the
                # direction, the phase and the evidence: none of those may collapse into one another
                "interaction": f"{e.direction} {e.phase} {e.evidence}",
                "direction": e.direction, "phase": e.phase, "evidence": e.evidence, "amount": e.amount,
                "change": e.change, "sd": e.sd, "n": e.n, "n_experiments": e.n_experiments, "p_value": e.p_value,
                "q_value": e.q_value,
                "window_start": e.window_start, "window_end": e.window_end, "exponential_h": e.exponential_h,
                "medium": e.medium,
                "study_ids": " ".join(e.study_ids), "experiments": " ".join(e.experiments),
                "cautions": " ".join(e.cautions), "notes": "; ".join(e.notes), "merged_arcs": e.merged_arcs,
                "merged_taxa": " ".join(e.merged_taxa), "line_style": line_style(e),
                "display_weight": display_weight(e)}
        # a null number becomes 0.0 in Cytoscape, which would read as a measured zero: leave it out
        edges.append({"data": {k: v for k, v in data.items() if v is not None}})
    return {"data": {"name": name, "shared_name": name}, "elements": {"nodes": nodes, "edges": edges}}


def _discrete(column: str, prop: str, mapping: dict, kind: str = "String") -> dict:
    return {"mappingType": "discrete", "mappingColumn": column, "mappingColumnType": kind,
            "visualProperty": prop,
            "map": [{"key": key, "value": value} for key, value in mapping.items()]}


def style(name: str = STYLE_NAME) -> dict:
    """The visual style, as CyREST takes it: defaults plus one mapping per visual property."""
    colors = {"produced": brand.PRODUCED, "consumed": brand.CONSUMED}
    return {
        "title": name,
        "defaults": [
            {"visualProperty": "NODE_SHAPE", "value": "ELLIPSE"},
            {"visualProperty": "NODE_FILL_COLOR", "value": brand.TAXON_NODE},
            {"visualProperty": "NODE_BORDER_PAINT", "value": brand.MUTED},
            {"visualProperty": "NODE_BORDER_WIDTH", "value": 1},
            {"visualProperty": "NODE_LABEL_COLOR", "value": brand.INK},
            {"visualProperty": "NODE_LABEL_FONT_SIZE", "value": 12},
            {"visualProperty": "NODE_LABEL_POSITION", "value": "S,N,c,0.00,4.00"},
            {"visualProperty": "NODE_SIZE", "value": 40},
            {"visualProperty": "EDGE_TRANSPARENCY", "value": 220},
            {"visualProperty": "EDGE_WIDTH", "value": 2},
            {"visualProperty": "EDGE_TARGET_ARROW_SHAPE", "value": "ARROW"},
        ],
        "mappings": [
            {"mappingType": "passthrough", "mappingColumn": "name", "mappingColumnType": "String",
             "visualProperty": "NODE_LABEL"},
            {"mappingType": "passthrough", "mappingColumn": "color", "mappingColumnType": "String",
             "visualProperty": "NODE_FILL_COLOR"},
            _discrete("kind", "NODE_SHAPE", {"taxon": "ELLIPSE", "metabolite": "ROUND_RECTANGLE"}),
            # metabolites smaller than taxa, so the taxa read as the hubs they are
            _discrete("kind", "NODE_SIZE", {"taxon": "40.0", "metabolite": "26.0"}),
            _discrete("direction", "EDGE_STROKE_UNSELECTED_PAINT", colors),
            _discrete("direction", "EDGE_TARGET_ARROW_UNSELECTED_PAINT", colors),
            _discrete("line_style", "EDGE_LINE_TYPE",
                      {"SOLID": "SOLID", "LONG_DASH": "LONG_DASH", "DOT": "DOT", "DASH_DOT": "DASH_DOT"}),
            # a whole-run change is no phase's: between the two, so it reads apart from both
            _discrete("phase", "EDGE_TRANSPARENCY", {"exponential": "230", "stationary": "110", "window": "230",
                                                     "whole_run": "170"}),
            {"mappingType": "continuous", "mappingColumn": "display_weight",
             "mappingColumnType": "Double", "visualProperty": "EDGE_WIDTH",
             "points": [{"value": 0.0, "lesser": "1.0", "equal": "1.0", "greater": "1.0"},
                        {"value": WIDE_AT_MM, "lesser": "12.0", "equal": "12.0", "greater": "12.0"}]},
        ],
    }


# Cytoscape's XML names for CyREST's column types
_XML_TYPES = {"String": "string", "Double": "float", "Integer": "integer", "Long": "long", "Boolean": "boolean"}


def style_xml(name: str = STYLE_NAME) -> str:
    """The same style as `style`, in the XML that File, Import, Styles from File reads.

    Cytoscape reads a style file only in this form: the CyREST JSON that `send` posts, and even the JSON
    Cytoscape itself exports, are refused with "Don't know how to read file" (Karoline, 3.10.4, 2026-09-28).
    """
    s = style(name)
    sections = {"NODE": [], "EDGE": []}
    props = {}
    for d in s["defaults"]:
        props[d["visualProperty"]] = ET.Element("visualProperty", name=d["visualProperty"],
                                                default=str(d["value"]))
    for m in s["mappings"]:
        vp = m["visualProperty"]
        prop = props.setdefault(vp, ET.Element("visualProperty", name=vp))
        kind = m["mappingType"]
        mapping = ET.SubElement(prop, f"{kind}Mapping", attributeName=m["mappingColumn"],
                                attributeType=_XML_TYPES[m["mappingColumnType"]])
        if kind == "discrete":
            for entry in m["map"]:
                ET.SubElement(mapping, "discreteMappingEntry", attributeValue=entry["key"], value=entry["value"])
        elif kind == "continuous":
            for p in m["points"]:
                ET.SubElement(mapping, "continuousMappingPoint", attrValue=str(p["value"]),
                              equalValue=p["equal"], greaterValue=p["greater"], lesserValue=p["lesser"])
    for vp, prop in props.items():
        sections[vp.split("_", 1)[0]].append(prop)
    vizmap = ET.Element("vizmap", id=f"VizMap-{name}", documentVersion="3.1")
    visual_style = ET.SubElement(vizmap, "visualStyle", name=name)
    ET.SubElement(visual_style, "network")
    node = ET.SubElement(visual_style, "node")
    ET.SubElement(node, "dependency", name="nodeSizeLocked", value="true")
    node.extend(sections["NODE"])
    edge = ET.SubElement(visual_style, "edge")
    # the arrowhead takes its color from its own mapping, the same as the line's
    ET.SubElement(edge, "dependency", name="arrowColorMatchesEdge", value="false")
    edge.extend(sections["EDGE"])
    ET.indent(vizmap, space="    ")
    return '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n' + ET.tostring(vizmap, "unicode") + "\n"


class CytoscapeError(RuntimeError):
    """Cytoscape could not be reached or refused the request, with what to do about it."""


# What can go wrong on the wire: urllib's own errors, and http.client's for a listener that answers with
# something that is not HTTP (a database, a dev server on the wrong port), which is not a URLError and
# escaped as a traceback before (found by Craig's agent, grownet #70). TimeoutError and OSError cover a socket
# that times out or resets outside urllib's wrapping.
WIRE_ERRORS = (urllib.error.URLError, http.client.HTTPException, TimeoutError, OSError)


def _local(url: str) -> str:
    """Every request goes to this machine: checked here, where the request is made, not only by callers."""
    if not url.startswith("http://127.0.0.1:"):
        raise CytoscapeError(f"refusing to send to {url!r}: foodnet only talks to a Cytoscape on this machine")
    return url


def _get(url: str, timeout: float = 30.0):
    with urllib.request.urlopen(_local(url), timeout=timeout) as response:     # noqa: S310 - checked
        return json.loads(response.read().decode("utf-8") or "null")


def _post(url: str, payload, timeout: float = 30.0):
    body = json.dumps(payload).encode("utf-8")
    request = urllib.request.Request(_local(url), data=body, method="POST",
                                     headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(request, timeout=timeout) as response:   # noqa: S310 - localhost only
        text = response.read().decode("utf-8")
    return json.loads(text) if text.strip() else {}


# Numbers an arc may not have: a presence_only arc has no amount, a single replicate no spread, a boolean
# search no numbers at all. They are left out of the payload rather than sent as null, because Cytoscape
# reads a null number as 0.0, which would read as a measured zero. The column is then declared here
# instead, so every arc carries all of them, with the cells of the arcs that have no value genuinely empty
# (grownet's finding, Karoline, 2026-10-03; checked against Cytoscape 3.10.3).
OPTIONAL_EDGE_COLUMNS = (("amount", "Double"), ("change", "Double"), ("sd", "Double"), ("n", "Integer"),
                         ("n_experiments", "Integer"),
                         ("p_value", "Double"), ("q_value", "Double"), ("window_start", "Double"),
                         ("window_end", "Double"), ("exponential_h", "Double"), ("merged_arcs", "Integer"))


def _declare_columns(root: str, suid, timeout: float) -> list:
    """Create the edge columns an arc may not carry, so the table holds every one. Returns their names."""
    table = f"{root}/networks/{suid}/tables/defaultedge"
    try:
        listed = _get(f"{table}/columns", timeout)
    except (urllib.error.URLError, http.client.HTTPException, ValueError):
        return []                      # the network is in; a missing column is not worth failing over
    if not isinstance(listed, list):
        listed = []                    # anything else is not a column list, so take none as present
    present = {column.get("name") for column in listed if isinstance(column, dict)}
    declared = []
    for column, kind in OPTIONAL_EDGE_COLUMNS:
        if column in present:
            continue
        try:
            _post(f"{table}/columns", {"name": column, "type": kind, "immutable": False}, timeout)
            declared.append(column)
        except (urllib.error.URLError, http.client.HTTPException, ValueError):
            continue
    return declared


def send(net: FoodNetwork, port: int = PORT, name: str = "foodnet",
         apply_style: bool = True, layout: str = "force-directed", timeout: float = 30.0) -> dict:
    """Post `net` into a running Cytoscape and return {"suid", "style", "warning", "columns", "url"}.

    Raises CytoscapeError with what to do when Cytoscape is not running, the port is wrong, or CyREST
    refuses the request. Nothing is sent anywhere but this port on the loopback interface.
    """
    if net.meta.get("incomplete"):
        name = f"{name} (INCOMPLETE)"        # records could not be read: the network says so where it lands
    root = base_url(port)
    try:
        created = _post(f"{root}/networks?title={urllib.parse.quote(name)}&collection={brand.NAME}",
                        network_json(net, name), timeout)
    except urllib.error.HTTPError as e:
        raise CytoscapeError(f"Cytoscape refused the network ({e.code} {e.reason}). Update Cytoscape to 3.8 or "
                             "later, whose CyREST takes networks this way.") from None
    except http.client.HTTPException as e:
        raise CytoscapeError(
            f"something is listening on port {port}, but it is not Cytoscape: it did not answer the way "
            f"Cytoscape's CyREST does ({type(e).__name__}). Check which port Cytoscape uses (its cyrest.port "
            "property, 1234 unless changed) and that no other program holds it.") from None
    except WIRE_ERRORS as e:
        raise CytoscapeError(_unreachable(port, getattr(e, "reason", e))) from None
    suid = created.get("networkSUID", created.get("data", {}).get("networkSUID"))
    if suid is None:
        raise CytoscapeError(f"Cytoscape accepted the network but reported no network id: {created!r}")

    declared = _declare_columns(root, suid, timeout)
    applied, warning = "", ""
    if apply_style:
        try:
            _ensure_style(root, timeout)
            # applying a style or a layout is a GET in CyREST; a POST is refused (405), which is how the
            # style went missing before (grownet #76)
            _get(f"{root}/apply/styles/{urllib.parse.quote(STYLE_NAME)}/{suid}", timeout)
            applied = STYLE_NAME
            if layout:
                _get(f"{root}/apply/layouts/{urllib.parse.quote(layout)}/{suid}", timeout)
        except urllib.error.HTTPError as e:
            warning = f"the network is in Cytoscape, but its style could not be applied ({e.code} {e.reason})"
        except WIRE_ERRORS as e:
            warning = f"the network is in Cytoscape, but its style could not be applied ({getattr(e, 'reason', e)})"
    return {"suid": suid, "style": applied, "warning": warning, "columns": declared,
            "url": f"{root}/networks/{suid}"}


def _unreachable(port: int, reason) -> str:
    """What to do when nothing answers on the CyREST port, worded for the page and the command line."""
    if "timed out" in str(reason).lower():
        return (f"Cytoscape did not answer on port {port} in time: it may still be starting. Wait until its "
                "window has fully opened, then try again.")
    return (f"Cytoscape is not running on this machine, or not listening on port {port} ({reason}). Start "
            "Cytoscape, wait until its window has fully opened, then try again. Cytoscape listens on port 1234 "
            "unless its cyrest.port property (Edit, Preferences, Properties) says otherwise; the command line "
            "takes another port with --cytoscape-port.")


def _ensure_style(root: str, timeout: float) -> None:
    """Make the foodnet style in Cytoscape the one this version draws.

    A new style is posted. One that exists is updated in place, so it always matches the legend: posting
    again would make Cytoscape keep the old one and add a renamed copy (foodnet_0), and leaving it alone
    would keep an older version's look. The defaults are replaced with PUT; the mappings are deleted one
    visual property at a time and posted again, since CyREST refuses deleting them all at once (405, found
    live on 3.10.3) and does not say whether a POST over an existing mapping replaces it.
    """
    wanted = style(STYLE_NAME)
    name = urllib.parse.quote(STYLE_NAME)
    if STYLE_NAME not in (_get(f"{root}/styles", timeout) or []):
        _post(f"{root}/styles", wanted, timeout)
        return
    _request("PUT", f"{root}/styles/{name}/defaults", wanted["defaults"], timeout)
    for mapping in _get(f"{root}/styles/{name}/mappings", timeout) or []:
        prop = urllib.parse.quote(mapping["visualProperty"])
        _request("DELETE", f"{root}/styles/{name}/mappings/{prop}", None, timeout)
    _post(f"{root}/styles/{name}/mappings", wanted["mappings"], timeout)


def _request(method: str, url: str, payload, timeout: float = 30.0):
    body = None if payload is None else json.dumps(payload).encode("utf-8")
    request = urllib.request.Request(_local(url), data=body, method=method,
                                     headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(request, timeout=timeout) as response:   # noqa: S310 - localhost only
        return response.read()
