"""The help names every setting, every field and every legend value (grownet's rule)."""
import html
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
    assert set(PHASES) == {"exponential", "stationary", "window", "whole_run"}


def test_the_help_renders_with_the_token_in_the_style_link():
    page = help_page.render_help("TOKEN", DEFAULTS, ("A", "B"))
    assert "/foodnet_style.xml?token=TOKEN" in page and "Every setting" in page


def test_the_texts_users_read_name_no_particular_figure():
    # Karoline, 2026-10-06: the advanced settings named the figure the tool was checked against, "but it's a
    # generic tool"
    from pathlib import Path

    from foodnet import gui
    root = Path(__file__).resolve().parents[1]
    texts = [gui.render_form("T"), help_page.render_help("T", DEFAULTS, ("A",)), help_page.render_about(),
             (root / "README.md").read_text(encoding="utf-8"), (root / "r" / "README.md").read_text(encoding="utf-8")]
    for text in texts:
        assert "figure 3" not in text.lower() and "paper" not in text.lower()


def test_the_help_warns_that_monocultures_do_not_predict_a_community():
    # Karoline, 2026-10-06: "update the help to clearly warn about the fact that what species produce/consume
    # alone is not necessarily predictive of what they do in a community"
    page = help_page.render_help("T", DEFAULTS, ("A",))
    assert 'id="alone"' in page and 'href="#alone"' in page
    assert "What a taxon produces and consumes on its own is not necessarily what it does in a" in page
    assert "not as a prediction" in page


def test_the_page_says_the_matrices_answer_a_mouseover(client):
    # Karoline, 2026-10-06: "add a small remark above the matrices that the hover exists"
    from conftest import run

    from foodnet import gui
    page = gui.render_result("T", run(client))
    remark = page.index("Hover over a cell to see its value")
    assert page.index("<h2 id=\"matrices\">Consumed and produced</h2>") < remark < page.index("<svg", remark)


def test_the_help_names_every_caution():
    from foodnet.model import CAUTIONS
    assert all(c in help_page.EDGE_FIELDS["cautions"] for c in CAUTIONS)


def test_the_help_explains_every_caution_under_its_tier():
    # Karoline, 2026-10-09: "in the help, add a section that explains the cautions"
    import html

    from foodnet import help as help_page
    from foodnet.legend import CAUTION_TEXT
    from foodnet.model import CAUTION_TIERS, CAUTIONS
    page = help_page.render_help("T", {}, ("A",))
    section = page[page.index('id="cautions"'):page.index('id="matrices"')]
    heads = [section.index(h) for h in ("Tier 1", "Tier 2", "Tier 3", "No value")]
    assert heads == sorted(heads)
    for caution in CAUTIONS:
        at = section.index("<dt><code>" + caution.replace("_", "_<wbr>") + "</code></dt>")
        tier = CAUTION_TIERS.get(caution)
        assert heads[(tier or 4) - 1] < at < ([*heads, len(section)])[tier or 4], caution
        assert html.escape(CAUTION_TEXT[caution][1:20], quote=True) in section, caution
    assert 'href="#cautions"' in page
    # every tier 1 caution says what to do (a review: a new one would ship without advice), and no_phase,
    # which is no caution but reaches cautions.csv, is explained
    advice = section[section.index("<ul>"):section.index("</ul>")]
    assert all(f"<code>{c}</code>" in advice for c in CAUTIONS if CAUTION_TIERS.get(c) == 1)
    assert "<code>no_phase</code>" in section


def test_no_legend_line_runs_past_the_picture():
    # a review measured six caution lines past the 900 px width
    import re
    for line in re.findall(r"<text [^>]*>([^<]*)</text>", legend.legend_svg()):
        assert len(html.unescape(line)) <= legend.LINE, line
