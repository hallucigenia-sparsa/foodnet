"""The JSON Schema of a foodnet network, and a validator that needs no third-party package.

The schema is built from the model's own fields and vocabularies, and the shipped file
(`schema/metabolite_network.schema.json`) must equal it: the guardrail gate checks both, and that an emitted
network validates (grownet's contract rule). `validate_document` checks what the schema says, by hand, since
the tool has no runtime dependencies.
"""
from __future__ import annotations

import json
import os
from dataclasses import fields

from .model import (
    CAUTIONS,
    DIRECTIONS,
    EVIDENCE,
    IDENTITIES,
    KINDS,
    KNOWN_SCHEMAS,
    PHASES,
    SCHEMA,
    Edge,
    Node,
    Study,
)

SCHEMA_FILE = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
                           "schema", "metabolite_network.schema.json")

_NUMBER_OR_NULL = {"type": ["number", "null"]}
_STRINGS = {"type": "array", "items": {"type": "string"}}


def _edge_properties() -> dict:
    props = {
        "source": {"type": "string"}, "target": {"type": "string"},
        "direction": {"enum": list(DIRECTIONS)}, "phase": {"enum": list(PHASES)}, "evidence": {"enum": list(EVIDENCE)},
        "amount": {"type": ["number", "null"], "minimum": 0}, "change": _NUMBER_OR_NULL, "sd": _NUMBER_OR_NULL,
        "n": {"type": ["integer", "null"]}, "n_experiments": {"type": ["integer", "null"]},
        "p_value": _NUMBER_OR_NULL, "q_value": _NUMBER_OR_NULL,
        "window_start": _NUMBER_OR_NULL, "window_end": _NUMBER_OR_NULL, "exponential_h": _NUMBER_OR_NULL,
        "medium": {"type": "string"},
        "study_ids": {"type": "array", "items": {"type": "string"}, "minItems": 1},
        "experiments": _STRINGS, "cautions": {"type": "array", "items": {"enum": list(CAUTIONS)}},
        "notes": _STRINGS, "merged_arcs": {"type": ["integer", "null"]}, "merged_taxa": _STRINGS,
    }
    assert set(props) == {f.name for f in fields(Edge)}, "schema and Edge disagree"
    return props


def _node_properties() -> dict:
    props = {"id": {"type": "string"}, "kind": {"enum": list(KINDS)}, "name": {"type": "string"},
             "identity": {"enum": ["", *IDENTITIES]}, "taxon_id": {"type": "string"},
             "species": {"type": "string"}, "chebi_id": {"type": "string"}}
    assert set(props) == {f.name for f in fields(Node)}, "schema and Node disagree"
    return props


SCHEMA_DOC = {
    "$schema": "https://json-schema.org/draft/2020-12/schema",
    "$id": SCHEMA,
    "title": "foodnet metabolite network",
    "description": "A bipartite network of taxa and metabolites: produced arcs run from a taxon to a metabolite, "
                   "consumed arcs from a metabolite to a taxon, each in one growth phase, with the studies behind it.",
    "type": "object",
    "required": ["schema", "nodes", "edges", "studies"],
    "properties": {
        "schema": {"enum": list(KNOWN_SCHEMAS)},
        "meta": {"type": "object"},
        "nodes": {"type": "array", "items": {"type": "object", "required": ["id", "kind"],
                                             "properties": _node_properties(), "additionalProperties": False}},
        "edges": {"type": "array", "items": {"type": "object",
                                             "required": ["source", "target", "direction", "phase", "evidence",
                                                          "study_ids"],
                                             "properties": _edge_properties(), "additionalProperties": False}},
        "studies": {"type": "array", "items": {"type": "object", "required": ["id"],
                                               "properties": {f.name: {"type": "string"} for f in fields(Study)},
                                               "additionalProperties": False}},
    },
}


def schema_json(indent: int = 2) -> str:
    return json.dumps(SCHEMA_DOC, indent=indent) + "\n"


def _type_ok(value, kind) -> bool:
    kinds = kind if isinstance(kind, list) else [kind]
    for k in kinds:
        if k == "null" and value is None:
            return True
        if k == "string" and isinstance(value, str):
            return True
        if k == "integer" and isinstance(value, int) and not isinstance(value, bool):
            return True
        if k == "number" and isinstance(value, (int, float)) and not isinstance(value, bool):
            return True
        if k == "array" and isinstance(value, list):
            return True
        if k == "object" and isinstance(value, dict):
            return True
    return False


def _check(value, rule: dict, where: str, problems: list) -> None:
    if "enum" in rule and value not in rule["enum"]:
        problems.append(f"{where}: {value!r} not one of {rule['enum']}")
        return
    if "type" in rule and not _type_ok(value, rule["type"]):
        problems.append(f"{where}: {value!r} is not {rule['type']}")
        return
    if isinstance(value, (int, float)) and "minimum" in rule and value < rule["minimum"]:
        problems.append(f"{where}: {value} is below {rule['minimum']}")
    if isinstance(value, list):
        if len(value) < rule.get("minItems", 0):
            problems.append(f"{where}: needs at least {rule['minItems']} item(s)")
        for i, item in enumerate(value):
            if "items" in rule:
                _check(item, rule["items"], f"{where}[{i}]", problems)


def _objects(items, spec: dict, label: str, problems: list) -> None:
    for i, obj in enumerate(items):
        where = f"{label}[{i}]"
        if not isinstance(obj, dict):
            problems.append(f"{where}: not an object")
            continue
        for key in spec["required"]:
            if key not in obj:
                problems.append(f"{where}: missing {key}")
        for key, value in obj.items():
            rule = spec["properties"].get(key)
            if rule is None:
                problems.append(f"{where}: unknown field {key}")
            else:
                _check(value, rule, f"{where}.{key}", problems)


def validate_document(doc) -> list:
    """Every way `doc` departs from the schema, plus the structural checks of the model; [] when valid."""
    if not isinstance(doc, dict):
        return ["the document is not a JSON object"]
    problems = []
    props = SCHEMA_DOC["properties"]
    for key in SCHEMA_DOC["required"]:
        if key not in doc:
            problems.append(f"missing {key}")
    if "schema" in doc and doc["schema"] not in KNOWN_SCHEMAS:
        problems.append(f"schema {doc['schema']!r} is not one foodnet reads ({', '.join(KNOWN_SCHEMAS)})")
    for part in ("nodes", "edges", "studies"):
        if isinstance(doc.get(part), list):
            _objects(doc[part], props[part]["items"], part, problems)
        elif part in doc:
            problems.append(f"{part}: not a list")
    if not problems:
        from .model import FoodNetwork
        problems += FoodNetwork.from_dict(doc).validate()
    return problems


def validate_file(path: str) -> list:
    with open(path, encoding="utf-8") as f:
        return validate_document(json.load(f))
