"""What makes two experiments' media the same medium.

mGrowthDB names a medium per compartment, but does not report its composition systematically: an added sugar,
a removed carbon source or a supplement usually lives only in the experiment's description or name ("WC plus
mucin beads", "RI_BH -Ac", "supplemented with 2g/L trehalose"). Matching on the medium name alone therefore
puts different media together (Karoline, 2026-10-06: "check descriptions that suggest something altered the
medium or atmosphere and treat it as another medium if you do. This default can be switched off in the
advanced settings").

A medium's identity (`identity`) is:

  * the medium names of its compartments, without case, punctuation or a parenthesized abbreviation (so
    "Wilkins-Chalgren Anaerobe Broth (WC)" and "Wilkins-Chalgren Anaerobe Broth" are one medium);
  * with the strict rule on, every alteration the description or the name states (`alterations`): something
    added ("plus", "supplemented with", "with initial", "+X", "X added", an amount such as "1 mM galactose")
    or taken away ("without", "with no", "no supplied", "-X", "X-free", "instead of X"), with its amount when
    one is stated, since 0.1% and 0.75% linoleic acid are two media;
  * with the strict rule on, the atmosphere when it is recorded (`atmosphere`): the gas composition of the
    compartments. An experiment that records none is not taken as different; it joins the recorded variant of
    its medium (`assign_atmospheres`).

Each pattern below was written against the descriptions mGrowthDB holds (all 559 experiments, 2026-10-06, read
again after a review found two phrasings and the amounts missed), and `tests/test_media.py` holds those phrasings
as cases. A phrasing nobody has used yet will not be recognized: the report lists every medium found, so a
reader can see what was told apart and exclude an experiment by hand.
"""
from __future__ import annotations

import re

# words that follow an alteration but do not name the compound
_TRAILING = re.compile(r"\b(?:beads?|added|supplied|for \d.*|as .*|in .*|at \d.*|mixed .*)$")
_STOP = {"wc", "medium", "media", "broth", "the", "a", "an"}
_COMPOUND = r"(?P<compound>[a-z0-9][a-z0-9 \-()/'.%\u00b5\u03bc]*?)"
# a compound ends at a sentence end (not a decimal point), a comma, or a word that starts what follows
_END = r"(?=\.(?!\d)|[;,]|\bfor\b|\bas\b|\bin\b|\bstarting\b|\busing\b|\bused\b|\binstead\b|$)"
_UNIT = r"(?:mm|mmol|g/l|g l-1|mg/l|mg/ml|%|um|\u00b5m|\u03bcm)"
_AMOUNT_TEXT = r"\d+(?:\.\d+)?\s*" + _UNIT + r"\s+(?:of\s+)?"
_AMOUNT = r"(?P<amount>" + _AMOUNT_TEXT + ")"

ADDED = [
    re.compile(r"\bplus\s+" + _AMOUNT + "?" + _COMPOUND + _END),
    re.compile(r"\bsupplemented with\s+" + _AMOUNT + "?" + _COMPOUND + _END),
    re.compile(r"\bwith (?:initial|added|additional)\s+" + _AMOUNT + "?" + _COMPOUND + _END),
    re.compile(r"\b(?:and|with)\s+" + _AMOUNT + _COMPOUND + _END),
    re.compile(r"\benriched with\s+" + _COMPOUND + _END),
    re.compile(r"\bspiked with\s+" + _COMPOUND + _END),
    # "0.1mg/L pantothenate added in the culture" (study SMGDB00000019)
    re.compile(r"(?:^|[\s(])" + _AMOUNT + _COMPOUND + r"\s+(?:added|supplied)\b"),
    # "WC + 10 mM acetate"
    re.compile(r"\+\s*" + _AMOUNT + _COMPOUND + _END),
    # "20 mM fructose instead of glucose": the fructose is added (the glucose taken away, below)
    re.compile(r"(?:^|[\s(])" + _AMOUNT + _COMPOUND + r"\s+instead of\b"),
]
REMOVED = [
    re.compile(r"\bwithout (?:initial |added |additional )?" + _COMPOUND + r"(?=\bplus\b|" + _END[3:]),
    re.compile(r"\bwith no (?:additional |added )?" + _COMPOUND + _END),
    # "No supplied pantothenate in culture" (study SMGDB00000019)
    re.compile(r"\bno (?:supplied|added|additional)\s+" + _COMPOUND + _END),
    re.compile(r"\blacking\s+" + _COMPOUND + _END),
    re.compile(r"\bdepleted of\s+" + _COMPOUND + _END),
    re.compile(r"\binstead of\s+" + _COMPOUND + _END),
    # "glucose-free WC"
    re.compile(r"\b(?P<compound>[a-z][a-z0-9\-]*?)-free\b"),
]
# in an experiment's name: "+Ac", "+GlcNAc", " -Ac" (a minus only after a space, so "WC-GlcNAc" is a name);
# letters only, so a code after it ("SE+Glu_1e4", a starting density) does not split one medium into many.
# A "+" also joins co-culture members in some names ("At+Ms"); foodnet reads monocultures, where it does not
NAME_ADDED = re.compile(r"\+\s*([A-Za-z]+)")
NAME_REMOVED = re.compile(r"(?:^|\s)[-\u2212]\s*([A-Za-z]+)")


def _amount(text: str) -> str:
    """An amount as one spelling: "1 mM ", "1mM" and "1.0 mM" are "1mm"; "5 g l-1" is "5g/l"; the micro sign
    and the Greek mu are one."""
    if not text:
        return ""
    text = re.sub(r"\s+of\s*$", "", text.strip())
    found = re.match(r"(\d+(?:\.\d+)?)\s*(.*)$", text)
    if not found:
        return re.sub(r"\s+", "", text)
    number = f"{float(found.group(1)):g}"
    unit = re.sub(r"\s+", "", found.group(2)).replace("gl-1", "g/l").replace("\u03bc", "\u00b5")
    return number + unit


def _clean(text: str) -> tuple:
    """(amount, compound) of one listed part: "1.5 µm tbhq antioxidant (dissolved in dmso)" is ("1.5µm",
    "tbhq antioxidant"); a remark in parentheses is not the compound."""
    text = re.sub(r"\([^)]*\)", " ", text).split("(")[0].strip()
    amount = ""
    found = re.match(r"^\s*(" + _AMOUNT_TEXT + ")", text)
    if found:
        amount, text = _amount(found.group(1)), text[found.end():]
    text = _TRAILING.sub("", " ".join(text.split()).strip(" ,"))
    words = [w for w in re.split(r"[\s,]+", text) if w and w not in _STOP]
    return amount, " ".join(words)


def _tokens(sign: str, match) -> set:
    found = set()
    groups = match.groupdict()
    amount = _amount(groups.get("amount"))
    for i, part in enumerate(re.split(r"\s+and\s+|,", match.group("compound"))):
        own, name = _clean(part)
        if name:
            # an amount is part of the medium: 0.1% and 0.75% linoleic acid are two media (study
            # SMGDB00000014), as are 0.1 and 3.0 mg/L pantothenate (study SMGDB00000019)
            dose = own or (amount if i == 0 else "")
            found.add(sign + (f"{dose} " if dose else "") + name)
    return found


def alterations(exp: dict) -> tuple:
    """The alterations an experiment's description and name state, as sorted tokens ("+mucin", "-glucose",
    "+1mm galactose")."""
    found = set()
    description = " ".join((exp.get("description") or "").casefold().split())
    for pattern in ADDED:
        for m in pattern.finditer(description):
            found |= _tokens("+", m)
    for pattern in REMOVED:
        for m in pattern.finditer(description):
            found |= _tokens("-", m)
    name = exp.get("name") or ""
    for m in NAME_ADDED.finditer(name):
        found.add("+" + _clean(m.group(1).casefold())[1])
    for m in NAME_REMOVED.finditer(name):
        found.add("-" + _clean(m.group(1).casefold())[1])
    return tuple(sorted(t for t in found if len(t) > 1))


def base_key(exp: dict) -> str:
    """The compartments' medium names, reduced to what tells media apart."""
    parts = []
    for c in exp.get("compartments", []):
        text = re.sub(r"\([^)]*\)", " ", (c.get("mediumName") or "").casefold())
        text = " ".join(re.sub(r"[^0-9a-z]+", " ", text).split())
        if text:
            parts.append(text)
    return "; ".join(sorted(set(parts))) or "unnamed medium"


GASES = ("O2", "CO2", "H2", "N2")


def atmosphere(exp: dict) -> str:
    """The recorded gas composition ("CO2 10, H2 10, N2 80"), or "" when no compartment records one."""
    parts = []
    for c in exp.get("compartments", []):
        values = [(g, c.get(g)) for g in GASES if c.get(g) not in (None, "")]
        if values:
            parts.append(", ".join(f"{g} {float(v):g}" for g, v in values if float(v) > 0) or "no gas recorded")
    return "; ".join(sorted(set(parts)))


def identity(exp: dict, strict: bool = True) -> dict:
    """{"key", "label", "alterations", "atmosphere"}: the medium an experiment ran in.

    `key` groups experiments into media; `label` is what a reader sees: the medium name with its alterations.
    With `strict` False only the medium names count, as before the strict rule."""
    from .reading import medium_of
    name = medium_of(exp) or "unnamed medium"
    if not strict:
        return {"key": base_key(exp), "label": name, "alterations": (), "atmosphere": ""}
    changed = alterations(exp)
    # spaces and hyphens inside a compound's name do not make another medium ("N-acetyl glucosamine")
    tokens = sorted({t[0] + re.sub(r"[\s\-]+", "", t[1:]) for t in changed})
    key = base_key(exp) + ("" if not tokens else " | " + " ".join(tokens))
    label = name + ("" if not changed else " (" + ", ".join(changed) + ")")
    return {"key": key, "label": label, "alterations": changed, "atmosphere": atmosphere(exp)}


def assign_atmospheres(identities: list) -> list:
    """The final keys and labels for a list of `identity` results: experiments of one medium that record
    different atmospheres become different media; one that records none joins the most common recorded variant
    of its medium (or stays as is when there is none)."""
    variants = {}
    for ident in identities:
        if ident["atmosphere"]:
            counts = variants.setdefault(ident["key"], {})
            counts[ident["atmosphere"]] = counts.get(ident["atmosphere"], 0) + 1
    out = []
    for ident in identities:
        counts = variants.get(ident["key"], {})
        if len(counts) <= 1:
            out.append(ident)          # one atmosphere, or none recorded: the medium is not split
            continue
        gas = ident["atmosphere"] or max(sorted(counts), key=lambda g: counts[g])
        out.append({**ident, "key": f"{ident['key']} | {gas}", "label": f"{ident['label']} [{gas}]"})
    return out
