"""Which mGrowthDB metabolite records are the same compound, and how their units become mM.

**One compound, two ChEBI ids.** mGrowthDB records an organic acid sometimes as the acid and sometimes as
its conjugate base: study SMGDB00000009 writes "acetic acid" (CHEBI:15366), study SMGDB00000007 "acetate"
(CHEBI:30089). At culture pH they are one pool, and an HPLC measures that pool whatever the record calls
it, so foodnet joins them under the conjugate base. Without this, the same compound measured in two studies
would be two metabolite nodes, and a taxon's acetate in one medium would never meet its acetate in another.
The table is explicit and short on purpose: a pair is joined only when it is listed here, with its reason,
and the node keeps every name and id it was recorded under.

(S)-lactic acid (CHEBI:422) joins lactate (CHEBI:24996): lactate in mGrowthDB is recorded without a
stereo-descriptor, and no study measures both forms side by side, so keeping them apart would only split one
measurement.

**Units.** mGrowthDB already converts most concentrations to mM (`techniqueUnits`), whatever they were
recorded in. A unit that is a concentration in another scale is converted here; g/L and mg/L only for a
compound whose molar mass is in MOLAR_MASS. Anything else (study SMGDB00000019 records "AUC", a peak area)
is not a concentration, and the series is left out and reported rather than read as one.
"""
from __future__ import annotations

# acid (or stereo-specific form) ChEBI id -> (the id it joins, why)
EQUIVALENT = {
    15366: (30089, "acetic acid is the conjugate acid of acetate"),
    30751: (15740, "formic acid is the conjugate acid of formate"),
    30772: (17968, "butyric acid is the conjugate acid of butyrate"),
    30768: (17272, "propionic acid is the conjugate acid of propionate"),
    16135: (48944, "isobutyric acid is the conjugate acid of isobutyrate"),
    # isovalerate is CHEBI:48942, as mGrowthDB records it (foodnet 0.2.0 joined it to 50128, biflavonoid, so
    # isovaleric acid and isovalerate stayed apart; found checking every id here against ChEBI, 2026-10-08)
    28484: (48942, "isovaleric acid is the conjugate acid of isovalerate"),
    422: (24996, "(S)-lactic acid is a form of lactate; mGrowthDB's lactate carries no stereo-descriptor"),
    15361: (15361, ""),        # pyruvate is recorded as the base already (listed so the name below applies)
    # ids checked against ChEBI (OLS, 2026-10-06), added after a review: an acid recorded beside its base would
    # otherwise be two metabolites (valeric acid and valerate are both in mGrowthDB, 2026-10-08)
    15741: (26806, "succinic acid is the conjugate acid of succinate"),
    30031: (26806, "succinate(2-) is the dianion of succinate"),
    17418: (31011, "valeric acid is the conjugate acid of valerate"),
}
# compounds that leave a culture as gas or vapor: their fall is no evidence of uptake (a thirteenth review round:
# ethanol evaporating while lactate concentrated passed as "used and made" in a culture that did not grow)
VOLATILE = {"ethanol", "methanol", "acetone", "acetaldehyde", "propanol", "1-propanol", "2-propanol",
            "isopropanol", "butanol", "1-butanol", "hydrogen", "dihydrogen", "carbon dioxide", "methane",
            "hydrogen sulfide", "dimethyl sulfide"}


# the same by ChEBI id, for a study that records one under another name: ethanol, methanol, acetone,
# acetaldehyde, 1-propanol, 2-propanol, 1-butanol, dihydrogen, carbon dioxide, methane, hydrogen sulfide
VOLATILE_CHEBI = {"16236", "17790", "15347", "15343", "28831", "17824", "28885", "18276", "16526", "16183", "16136"}


def volatile(name: str, chebi_id: str = "") -> bool:
    return (" ".join((name or "").casefold().split()) in VOLATILE
            or str(chebi_id or "").removeprefix("CHEBI:") in VOLATILE_CHEBI)


# the name a joined compound is shown under
NAMES = {30089: "acetate", 15740: "formate", 17968: "butyrate", 17272: "propionate", 48944: "isobutyrate",
         48942: "isovalerate", 24996: "lactate", 15361: "pyruvate", 26806: "succinate", 31011: "valerate"}

# concentration units to mM
SCALE_TO_MM = {"mm": 1.0, "mmol/l": 1.0, "mmol l-1": 1.0, "mmol/litre": 1.0,
               "um": 1e-3, "µm": 1e-3, "μm": 1e-3, "umol/l": 1e-3, "µmol/l": 1e-3, "μmol/l": 1e-3,
               "nm": 1e-6, "nmol/l": 1e-6, "m": 1000.0, "mol/l": 1000.0}
# mass concentration units to g/L
MASS_TO_G_PER_L = {"g/l": 1.0, "mg/l": 1e-3, "mg/ml": 1.0, "ug/ml": 1e-3, "µg/ml": 1e-3}
# g/mol, by the ChEBI id foodnet keys the compound under (after EQUIVALENT)
MOLAR_MASS = {17234: 180.156, 15361: 88.06, 27082: 342.296, 30089: 60.052, 15740: 46.025, 17968: 88.106,
              17272: 74.079, 48944: 88.106, 48942: 102.133, 24996: 90.078, 26806: 118.088, 28757: 180.156,
              28260: 180.156, 37684: 180.156, 18222: 150.13, 17057: 342.297, 33984: 164.16, 31011: 102.133}


def canonical(chebi_id, name: str) -> dict:
    """{"id", "chebi_id", "name", "joined"}: the node a metabolite record belongs to.

    `id` is "chebi:<id>" after joining equivalent forms, or "metabolite:<name>" when mGrowthDB gives no
    ChEBI id. `joined` is the reason the record's own id was replaced, or ""."""
    try:
        # ChEBI writes its ids with a "CHEBI:" prefix, which mGrowthDB leaves off
        chebi = (int(str(chebi_id).strip().upper().removeprefix("CHEBI:"))
                 if chebi_id not in (None, "") else None)
    except (TypeError, ValueError):
        chebi = None
    if chebi is None:
        key = " ".join((name or "unnamed").casefold().split())
        return {"id": f"metabolite:{key}", "chebi_id": "", "name": name or key, "joined": ""}
    target, why = EQUIVALENT.get(chebi, (chebi, ""))
    return {"id": f"chebi:{target}", "chebi_id": str(target), "name": NAMES.get(target, name),
            "joined": why if target != chebi else ""}


def to_mm(unit: str, chebi_id: str = "") -> tuple:
    """(factor that turns a value in `unit` into mM, None) or (None, why it cannot be converted)."""
    key = " ".join((unit or "").strip().casefold().split())
    if key in SCALE_TO_MM:
        return SCALE_TO_MM[key], None
    if key in MASS_TO_G_PER_L:
        try:
            mass = MOLAR_MASS[int(chebi_id)]
        except (KeyError, TypeError, ValueError):
            return None, f"recorded in {unit}, and foodnet has no molar mass for this compound to convert it"
        return MASS_TO_G_PER_L[key] / mass * 1000.0, None
    if not key:
        return None, "no unit recorded"
    return None, f"recorded in {unit!r}, which is not a concentration foodnet can convert to mM"
