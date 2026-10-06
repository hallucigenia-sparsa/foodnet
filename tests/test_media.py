"""What makes two experiments' media one medium (foodnet.media).

Karoline, 2026-10-06: "check descriptions that suggest something altered the medium or atmosphere and treat it
as another medium if you do. This default can be switched off in the advanced settings." The descriptions below
are phrasings mGrowthDB holds (2026-10-06), so the rules are checked on what they will meet.
"""
import pytest

from foodnet import media


def _exp(description, name="x", medium="Wilkins-Chalgren Anaerobe Broth (WC)", **gas):
    return {"description": description, "name": name, "compartments": [{"mediumName": medium, **gas}]}


@pytest.mark.parametrize("description, name, tokens", [
    ("BT with WC plus mucin beads for 120 h", "BT_MUCIN", ("+mucin",)),
    ("BT with WC for 120h", "BT_WC", ()),
    ("RI and BH co-culture with initial acetate", "RI_BH +Ac", ("+ac", "+acetate")),
    ("RI and BH co-culture without initial acetate", "RI_BH -Ac", ("-ac", "-acetate")),
    ("WC plus 1 mM galactose", "x", ("+1mm galactose",)),
    ("WC and 1mM N-acetylglucosamine", "x", ("+1mm n-acetylglucosamine",)),
    ("modified WC without glucose and pyruvate plus 1 mM mannose", "x", ("+1mm mannose", "-glucose", "-pyruvate")),
    ("C. difficile R20291 (DSM 27147) monoculture grown in DM29 supplemented with 2g/L trehalose.", "x",
     ("+2g/l trehalose",)),
    ("C. difficile R20291 (DSM 27147) monoculture grown in DM29 with no additional carbohydrates added.", "x",
     ("-carbohydrates",)),
    ("monoculture grown in modified 10% BHI supplemented with Glutamate. Starting density of ~10^4 CFUs/mL.",
     "SA+Glu_1e4", ("+glu", "+glutamate")),
    ("At monoculture grown on a minimal medium with 0.75% linoleic acid, used for measuring reactive oxygen", "x",
     ("+0.75% linoleic acid",)),
    ("At monoculture grown on a minimal medium with 0.75% linoleic acid and 1.5μm TBHQ antioxidant "
     "(dissolved in DMSO)", "x", ("+0.75% linoleic acid", "+1.5μm tbhq antioxidant")),
    # found missed by a review (2026-10-06): an amount "added", "no supplied", and phrasings to expect
    ("0.1mg/L pantothenate added in the culture", "x", ("+0.1mg/l pantothenate",)),
    ("No supplied pantothenate in culture", "x", ("-pantothenate",)),
    ("WC + 10 mM acetate", "x", ("+10mm acetate",)),
    ("glucose-free WC", "x", ("-glucose",)),
    ("WC with 20 mM fructose instead of glucose", "x", ("+20mm fructose", "-glucose")),
    ("E. coli LF82 grown in WC as monoculture for 168h", "EC", ()),
    ('Roseburia monoculture controls of the "btri" experiment', "ri2", ()),
])
def test_alterations_stated_in_descriptions_and_names(description, name, tokens):
    assert media.alterations(_exp(description, name)) == tokens


def test_amounts_tell_media_apart():
    low = media.identity(_exp("minimal medium with 0.1% linoleic acid"))
    high = media.identity(_exp("minimal medium with 0.75% linoleic acid"))
    same = media.identity(_exp("minimal medium with 0.1 % linoleic acid"))
    assert low["key"] != high["key"] and low["key"] == same["key"]


def test_spellings_of_one_compound_are_one_medium():
    a = media.identity(_exp("modified WC without glucose and pyruvate plus 1 mM N-acetyl glucosamine"))
    b = media.identity(_exp("modified WC without glucose and pyruvate plus 1 mM N-acetylglucosamine"))
    assert a["key"] == b["key"] and a["key"] != media.identity(_exp("WC plus 1 mM N-acetylglucosamine"))["key"]


def test_the_strict_rule_can_be_switched_off():
    plain, mucin = _exp("BT with WC for 120h"), _exp("BT with WC plus mucin beads for 120 h")
    assert media.identity(plain)["key"] != media.identity(mucin)["key"]
    assert media.identity(plain, strict=False)["key"] == media.identity(mucin, strict=False)["key"]


def test_a_recorded_atmosphere_splits_a_medium_and_an_unrecorded_one_does_not():
    a = media.identity(_exp("x", CO2="10.00", H2="10.00", N2="80.00"))
    b = media.identity(_exp("x", CO2="20.00", N2="80.00"))
    unknown = media.identity(_exp("x"))
    keys = [i["key"] for i in media.assign_atmospheres([a, a, b, unknown])]
    assert keys[0] == keys[1] != keys[2]
    assert keys[3] == keys[0]                     # joins the most common recorded variant
    # one recorded atmosphere and some unrecorded ones: one medium, not split
    assert len({i["key"] for i in media.assign_atmospheres([a, unknown, unknown])}) == 1
