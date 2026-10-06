"""The foodnet network: which taxon produces and which consumes which metabolite.

A FoodNetwork is bipartite. Its nodes are taxa (strains, as mGrowthDB records them) and metabolites (by
ChEBI id), and every arc joins one of each:

  * a **produced** arc runs from the taxon to the metabolite: the compound rose in the taxon's monoculture;
  * a **consumed** arc runs from the metabolite to the taxon: the compound fell.

Each arc belongs to one growth phase (exponential, stationary, or a time window the user set), since what a
culture makes or takes up while it grows can differ from what it does once it has stopped growing.

An arc is either **measured**, with an amount (the net change in mM over the phase, averaged over
replicates), or **presence_only**: the change was seen in another medium than the one the values come from,
so its direction counts and its size does not (Karoline, 2026-10-04). A compound that was never assayed has
no arc and no number; it is never written as zero.

Every arc carries the studies it rests on, so attribution resolves at the arc level, as in grownet.
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field

SCHEMA = "foodnet.metabolite_network/v0"
KNOWN_SCHEMAS = (SCHEMA,)

KINDS = ("taxon", "metabolite")
DIRECTIONS = ("produced", "consumed")
PHASES = ("exponential", "stationary", "window")
EVIDENCE = ("measured", "presence_only")
# what a node's id rests on: the NCBI taxon id of the strain, genus and species of its name, a genus (after
# merging to genus), or the ChEBI id of a metabolite (or its name, when mGrowthDB gives no ChEBI id)
IDENTITIES = ("ncbi", "name", "genus", "chebi", "metabolite_name")
# what a reader should know about an arc, without it being wrong
CAUTIONS = (
    "single_replicate",        # one replicate: no spread
    "short_record",            # the metabolite series covers less than 24 h
    "window_beyond_data",      # the phase or window starts before the first or ends after the last sample
    "stationary_not_reached",  # the growth curve was still rising at its last point: no stationary phase
    "conflict",                # experiments (or, merged to genus, taxa) disagree: no value, no arc
    "not_detected_in_value_medium",  # seen in another medium, assayed and not seen in the value medium
    "phase_from_other_replicates",   # this replicate had no growth curve; the boundary of its siblings is used
    "coarse_sampling",         # the phase boundary rests on fewer than three growth samples
    "boundaries_differ",       # the replicates' phase boundaries lie further apart than a sampling interval
    "no_variance",             # the replicates are identical (rounding, or one series deposited twice): no test
    "amounts_differ",          # the experiments agree in direction, with amounts more than twofold apart
    "experiment_left_out",     # an inconclusive experiment that does not contradict the others is left out
    "still_changing",          # the compound kept changing beyond the limit right after growth slowed
    "growth_rate_boundary",    # the growth-rate rule ended growth more than one sample before the 90% rule
    "start_differs",           # the replicates' windows start more than a quarter of the phase apart
    "pair_decided",            # an experiment of two replicates decided it (a change or no change), no interval
)

# words before a genus that are not a genus (the same rule as grownet)
GENUS_QUALIFIERS = ("candidatus", "unclassified", "uncultured")


def genus_name(name: str) -> str:
    """The genus of an organism name, as mGrowthDB writes it: its first word after any qualifier. NCBI's
    brackets stay, so "[Clostridium] scindens" is placed apart from Clostridium."""
    words = (name or "").split()
    while words and words[0].lower() in GENUS_QUALIFIERS:
        words = words[1:]
    if not words:
        return "unknown"
    first = words[0]
    if first.startswith("[") and first.endswith("]"):
        return "[" + first[1:-1].capitalize() + "]"
    return first.capitalize()


def genus_species(name: str) -> str:
    """Genus and species key for matching a strain across experiments (drops the strain designation)."""
    return " ".join((name or "").split()[:2]).lower()


@dataclass(frozen=True)
class Study:
    """A source study behind one or more arcs (the unit of attribution)."""

    id: str
    citation: str = ""
    license: str = ""
    url: str = ""


@dataclass(frozen=True)
class Node:
    """A taxon or a metabolite."""

    id: str
    kind: str                     # one of KINDS
    name: str = ""
    identity: str = ""            # one of IDENTITIES
    taxon_id: str = ""            # a taxon's NCBI taxon id, as mGrowthDB records it
    species: str = ""             # a taxon's genus and species, from its name
    chebi_id: str = ""            # a metabolite's ChEBI id


@dataclass(frozen=True)
class Edge:
    """A produced (taxon to metabolite) or consumed (metabolite to taxon) arc in one growth phase."""

    source: str
    target: str
    direction: str                # one of DIRECTIONS
    phase: str                    # one of PHASES
    evidence: str                 # one of EVIDENCE
    amount: float | None = None   # mM, the size of the net change in this direction; None for presence_only
    change: float | None = None   # mM, the signed mean net change (positive = produced)
    sd: float | None = None       # mM, standard deviation over the experiments' means (over replicates in one)
    n: int | None = None          # replicates behind the mean
    n_experiments: int | None = None  # experiments behind the mean: the unit the mean and the test count
    p_value: float | None = None  # one-sample t-test of the experiment means (one: replicates) against zero
    q_value: float | None = None  # p_value corrected for multiple testing over the arcs' family of tests
    window_start: float | None = None  # h, mean start of the phase over the replicates
    window_end: float | None = None    # h, mean end of the phase over the replicates
    exponential_h: float | None = None  # h, how long the cultures grew exponentially (mean over replicates)
    medium: str = ""              # the medium (or media) the arc rests on
    study_ids: tuple = ()
    experiments: tuple = ()
    cautions: tuple = ()          # CAUTIONS
    notes: tuple = ()
    merged_arcs: int | None = None      # arcs merged into this one (merge across studies), None when not merged
    merged_taxa: tuple = ()             # with merging to genus: the taxa behind the arc

    @property
    def taxon(self) -> str:
        return self.source if self.direction == "produced" else self.target

    @property
    def metabolite(self) -> str:
        return self.target if self.direction == "produced" else self.source

    def validate(self) -> list:
        problems = []
        where = f"arc {self.source}->{self.target}"
        if self.direction not in DIRECTIONS:
            problems.append(f"{where}: direction {self.direction!r} not in {DIRECTIONS}")
        if self.phase not in PHASES:
            problems.append(f"{where}: phase {self.phase!r} not in {PHASES}")
        if self.evidence not in EVIDENCE:
            problems.append(f"{where}: evidence {self.evidence!r} not in {EVIDENCE}")
        if self.evidence == "presence_only" and self.amount is not None:
            problems.append(f"{where}: a presence_only arc carries no amount")
        if self.amount is not None and self.amount < 0:
            problems.append(f"{where}: amount {self.amount} is negative")
        for flag in self.cautions:
            if flag not in CAUTIONS:
                problems.append(f"{where}: caution {flag!r} not in {CAUTIONS}")
        if not self.study_ids:
            problems.append(f"{where}: no study_ids (arc-level attribution requires at least one)")
        return problems


@dataclass
class FoodNetwork:
    nodes: dict = field(default_factory=dict)     # id -> Node
    edges: list = field(default_factory=list)     # list[Edge]
    studies: dict = field(default_factory=dict)   # id -> Study
    meta: dict = field(default_factory=dict)
    schema: str = SCHEMA

    def add_study(self, study: Study) -> None:
        self.studies[study.id] = study

    def add_node(self, node: Node) -> None:
        self.nodes[node.id] = node

    def add_edge(self, edge: Edge) -> None:
        self.edges.append(edge)

    def taxa(self) -> list:
        return [n for n in self.nodes.values() if n.kind == "taxon"]

    def metabolites(self) -> list:
        return [n for n in self.nodes.values() if n.kind == "metabolite"]

    def validate(self) -> list:
        problems = []
        for n in self.nodes.values():
            if n.kind not in KINDS:
                problems.append(f"node {n.id}: kind {n.kind!r} not in {KINDS}")
        for e in self.edges:
            problems += e.validate()
            for end in (e.source, e.target):
                if end not in self.nodes:
                    problems.append(f"arc references missing node {end!r}")
            if e.taxon in self.nodes and self.nodes[e.taxon].kind != "taxon":
                problems.append(f"arc {e.source}->{e.target}: its taxon end {e.taxon!r} is not a taxon")
            if e.metabolite in self.nodes and self.nodes[e.metabolite].kind != "metabolite":
                problems.append(f"arc {e.source}->{e.target}: its metabolite end {e.metabolite!r} is not a metabolite")
            for sid in e.study_ids:
                if sid not in self.studies:
                    problems.append(f"arc {e.source}->{e.target} cites unknown study {sid!r}")
        return problems

    def to_dict(self) -> dict:
        return {"schema": self.schema, "meta": self.meta,
                "nodes": [asdict(n) for n in self.nodes.values()],
                "edges": [asdict(e) for e in self.edges],
                "studies": [asdict(s) for s in self.studies.values()]}

    def to_json(self, indent: int = 2) -> str:
        return json.dumps(self.to_dict(), indent=indent, ensure_ascii=False)

    @classmethod
    def from_dict(cls, d: dict) -> FoodNetwork:
        net = cls(meta=d.get("meta", {}), schema=d.get("schema", SCHEMA))
        for s in d.get("studies", []):
            net.add_study(Study(**s))
        for n in d.get("nodes", []):
            net.add_node(Node(**n))
        for e in d.get("edges", []):
            e = dict(e)
            for key in ("study_ids", "experiments", "cautions", "notes", "merged_taxa"):
                e[key] = tuple(e.get(key, ()))
            net.add_edge(Edge(**e))
        return net
