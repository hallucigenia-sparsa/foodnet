"""A small synthetic mGrowthDB, and the network it gives, for every test that needs data.

Nothing here is real data (AGENTS.md: only synthetic fixtures belong in the repository). The world is
chosen so every number can be checked by hand:

Study SMGDB00000001, medium "Wilkins-Chalgren Anaerobe Broth (WC)":
  E1, Alpha alpha A1 (taxon 1), two replicates. Growth (culture flow cytometry) at 0, 4, 8, 12, 24 h:
      1, 10, 100, 1000, 1000. Start 1, maximum 1000, so exponential growth ends at the first sample at or
      above 1 + 0.9 * 999 = 900.1, which is 12 h.
      replicate 1: glucose 10, 8, 2, 2 and acetic acid 0, 1, 4, 5 at 0, 6, 12, 24 h
      replicate 2: glucose 10, 7, 2, 1 and acetic acid 0, 1, 4, 3
      exponential (0 to 12 h): glucose -8 and -8 (mean -8, consumed 8); acetate +4 and +4 (produced 4)
      stationary (12 to 24 h): glucose 0 and -1 (mean -0.5, consumed 0.5); acetate +1 and -1 (mean 0)
  E2, Beta beta B1 (taxon 2), one replicate. Growth (OD) at 0, 8, 16, 24 h: 0.01, 0.1, 0.5, 0.5, so the
      boundary is the first sample at or above 0.01 + 0.9 * 0.49 = 0.451: 16 h.
      glucose 10, 9, 4, 4 and butyric acid 0, 0, 3, 3 at 0, 8, 16, 24 h
      exponential: glucose -6 (consumed 6), butyrate +3 (produced 3); stationary: both 0

Study SMGDB00000002, medium "mMCB":
  E3, Beta beta B1, one replicate: growth 0.01, 0.1, 0.5, 0.5 at 0, 8, 16, 24 h; formate 0, 2, 5, 5
      (exponential +5: produced, but only in another medium than the values).
  E4, Gamma gamma C1 (taxon 3), one replicate: same growth; glucose 10, 8, 5, 5 (exponential -5).

Media: WC holds two taxa (A, B) and mMCB two (B, C), a tie. WC has three replicates against mMCB's two, so
WC gives the values, and the tie is reported.
"""
import copy

import pytest

STUDIES = {
    "SMGDB00000001": {"id": "SMGDB00000001", "name": "Synthetic study one", "url": "",
                      "uploadedAt": "2026-01-01", "publishedAt": "2026-01-02",
                      "experiments": [{"id": "EMGDB000000001", "name": "A_WC"},
                                      {"id": "EMGDB000000002", "name": "B_WC"}]},
    "SMGDB00000002": {"id": "SMGDB00000002", "name": "Synthetic study two", "url": "",
                      "uploadedAt": "2026-01-01", "publishedAt": "2026-01-03",
                      "experiments": [{"id": "EMGDB000000003", "name": "B_MCB"},
                                      {"id": "EMGDB000000004", "name": "C_MCB"}]},
}
STRAINS = {1: "Alpha alpha A1", 2: "Beta beta B1", 3: "Gamma gamma C1"}
CHEBI = {"glucose": 17234, "acetic acid": 15366, "butyric acid": 30772, "formate": 15740}


def _experiment(eid, study, taxon, medium, replicates):
    return {"id": eid, "name": eid, "description": f"{STRAINS[taxon]} in {medium}", "studyId": study,
            "cultivationMode": "batch",
            "communityStrains": [{"id": taxon, "NCBId": taxon, "name": STRAINS[taxon]}],
            "compartments": [{"name": "c", "mediumName": medium}],
            "bioreplicates": [{"id": rid, "name": rid} for rid in replicates]}


EXPERIMENTS = {
    "EMGDB000000001": _experiment("EMGDB000000001", "SMGDB00000001", 1, "Wilkins-Chalgren Anaerobe Broth (WC)",
                                  ["a1", "a2"]),
    "EMGDB000000002": _experiment("EMGDB000000002", "SMGDB00000001", 2, "Wilkins-Chalgren Anaerobe Broth (WC)",
                                  ["b1"]),
    "EMGDB000000003": _experiment("EMGDB000000003", "SMGDB00000002", 2, "mMCB", ["b2"]),
    "EMGDB000000004": _experiment("EMGDB000000004", "SMGDB00000002", 3, "mMCB", ["c1"]),
}

FAST = [(0, 1), (4, 10), (8, 100), (12, 1000), (24, 1000)]
SLOW = [(0, 0.01), (8, 0.1), (16, 0.5), (24, 0.5)]
SERIES = {
    "a1": {"growth": ("fc", FAST),
           "glucose": [(0, 10), (6, 8), (12, 2), (24, 2)], "acetic acid": [(0, 0), (6, 1), (12, 4), (24, 5)]},
    "a2": {"growth": ("fc", FAST),
           "glucose": [(0, 10), (6, 7), (12, 2), (24, 1)], "acetic acid": [(0, 0), (6, 1), (12, 4), (24, 3)]},
    "b1": {"growth": ("od", SLOW), "glucose": [(0, 10), (8, 9), (16, 4), (24, 4)],
           "butyric acid": [(0, 0), (8, 0), (16, 3), (24, 3)]},
    "b2": {"growth": ("od", SLOW), "formate": [(0, 0), (8, 2), (16, 5), (24, 5)]},
    "c1": {"growth": ("od", SLOW), "glucose": [(0, 10), (8, 8), (16, 5), (24, 5)]},
}


class FakeClient:
    """The reads foodnet makes, answered from the synthetic world above."""

    def __init__(self, series=None):
        self.series = copy.deepcopy(series or SERIES)
        self.contexts = {}
        self.bioreplicates = {}
        cid = 100
        for rid, data in self.series.items():
            contexts = []
            for name, value in data.items():
                cid += 1
                if name == "growth":
                    technique, points = value
                    contexts.append({"id": cid, "techniqueType": technique, "techniqueUnits": "",
                                     "subject": {"type": "bioreplicate", "name": rid}})
                else:
                    points = value
                    contexts.append({"id": cid, "techniqueType": "metabolite", "techniqueUnits": "mM",
                                     "subject": {"type": "metabolite", "name": name, "chebiId": CHEBI[name]}})
                self.contexts[cid] = points
            self.bioreplicates[rid] = {"id": rid, "name": rid, "isAverage": False, "measurementTimeUnits": "h",
                                       "measurementContexts": contexts}

    def get_study(self, sid):
        from foodnet.mgrowthdb import MGrowthDBError
        if sid not in STUDIES:
            raise MGrowthDBError(f"no study {sid}", status=404)
        return STUDIES[sid]

    def get_experiment(self, eid):
        return EXPERIMENTS[eid]

    def get_bioreplicate(self, rid):
        return self.bioreplicates[rid]

    def get_measurement_series(self, cid):
        return [(t, v, None) for t, v in self.contexts[cid]]

    def study_experiments(self, sid):
        return [EXPERIMENTS[e["id"]] for e in self.get_study(sid)["experiments"]]

    def search(self, **_):
        return {"studies": list(STUDIES)}


@pytest.fixture
def client():
    return FakeClient()


@pytest.fixture
def result(client):
    """The exponential-phase result for all three taxa, with growth rates."""
    from foodnet.search import run_query
    return run_query(client, ["Alpha alpha", "Beta beta", "Gamma gamma"], {"report_rates": True})


def run(client, **settings):
    from foodnet.search import run_query
    return run_query(client, ["Alpha alpha", "Beta beta", "Gamma gamma"], settings)
