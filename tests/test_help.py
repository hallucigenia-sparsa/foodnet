"""The help names every setting, every field and every legend value (grownet's rule)."""
from dataclasses import fields

from foodnet import help as help_page
from foodnet import legend
from foodnet.model import CAUTIONS, DIRECTIONS, EVIDENCE, PHASES, Edge, Node
from foodnet.search import DEFAULTS


def test_every_setting_has_a_help_line():
    assert set(help_page.SETTINGS) == set(DEFAULTS)


def test_every_field_has_a_help_line():
    assert set(help_page.EDGE_FIELDS) == {f.name for f in fields(Edge)}
    assert set(help_page.NODE_FIELDS) == {f.name for f in fields(Node)}


def test_the_legend_names_every_value():
    svg = legend.legend_svg()
    for word in (*DIRECTIONS, *EVIDENCE, *CAUTIONS, "stationary", "exponential", "window"):
        assert word in svg, word
    assert set(legend.CAUTION_TEXT) == set(CAUTIONS)
    assert set(PHASES) == {"exponential", "stationary", "window"}


def test_the_help_renders_with_the_token_in_the_style_link():
    page = help_page.render_help("TOKEN", DEFAULTS, ("A", "B"))
    assert "/foodnet_style.xml?token=TOKEN" in page and "Every setting" in page
