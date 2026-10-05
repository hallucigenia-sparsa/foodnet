"""The image of the two matrices (Karoline, 2026-10-04: "as the first result, generate and show a downloadable
image of the 2 matrices (as in Figure 3c)")."""
import xml.etree.ElementTree as ET

from conftest import run

from foodnet import figure

SVG = "{http://www.w3.org/2000/svg}"


def _texts(svg):
    return ["".join(t.itertext()) for t in ET.fromstring(svg).iter(f"{SVG}text")]


def test_both_panels_with_the_numbers_and_the_key(client):
    svg = figure.matrices_svg(run(client))
    texts = _texts(svg)
    assert "Consumed (mM)" in texts and "Produced (mM)" in texts
    assert "8.0" in texts and "4.0" in texts and "6.0" in texts and "3.0" in texts   # A, B as in conftest.py
    assert "not assayed" in texts and "seen only in another medium" in texts
    assert "Alpha alpha A1" in texts


def test_cell_states_have_their_own_marks(client):
    root = ET.fromstring(figure.matrices_svg(run(client)))
    fills = [r.get("fill") for r in root.iter(f"{SVG}rect")]
    assert figure.NOT_ASSAYED in fills                       # C's acetate, never assayed in WC
    # presence only: B's formate and C's glucose, one circle each, plus the one in the key
    assert len(list(root.iter(f"{SVG}circle"))) == 3


def test_booleans_print_no_numbers(client):
    texts = _texts(figure.matrices_svg(run(client, booleans=True)))
    assert "Consumed" in texts and not any(t.replace(".", "").isdigit() for t in texts)


def test_the_page_shows_it_first_and_offers_it_as_a_download(client):
    from foodnet import gui
    page = gui.render_result("T", run(client))
    assert page.index("Consumed and produced") < page.index("<h2>Taxa</h2>")
    assert "/matrices.svg?token=T" in page and "Download the image (.svg)" in page


def test_a_wide_matrix_stacks_its_panels(client):
    from foodnet import figure as f
    root = ET.fromstring(f.matrices_svg(run(client, phase="both")))
    titles = {"".join(t.itertext()): float(t.get("y")) for t in root.iter(f"{SVG}text")
              if "".join(t.itertext()) in ("Consumed (mM)", "Produced (mM)")}
    # 4 metabolites x 2 phases = 8 columns per panel: 352 px each, side by side about 840 px, under 1000
    assert titles["Consumed (mM)"] == titles["Produced (mM)"]
    old = f.MAX_SIDE_BY_SIDE
    f.MAX_SIDE_BY_SIDE = 500
    try:
        root = ET.fromstring(f.matrices_svg(run(client, phase="both")))
    finally:
        f.MAX_SIDE_BY_SIDE = old
    titles = {"".join(t.itertext()): float(t.get("y")) for t in root.iter(f"{SVG}text")
              if "".join(t.itertext()) in ("Consumed (mM)", "Produced (mM)")}
    assert titles["Produced (mM)"] > titles["Consumed (mM)"]


def test_the_name_in_the_caption_stays_text_inside_the_svg(client):
    from foodnet import gui
    page = gui.render_result("T", run(client))
    svg = page[page.index("<svg"):page.index("</svg>")]
    assert '<span class="name">' not in svg


def test_the_caption_wraps_inside_the_image():
    from foodnet import figure as f
    lines = f._wrap("one two three four five six seven", 60, 11)
    assert len(lines) > 1 and all(f._text_width(x, 11) <= 60 or " " not in x for x in lines)


def test_second_window_columns_carry_a_star_that_the_caption_explains(client):
    texts = _texts(figure.matrices_svg(run(client, second_window_metabolites="glucose")))
    assert "glucose *" in texts and "acetate" in texts          # no suffix on the columns of the phase
    assert any("(* 0 h to the last sample)" in t for t in texts)
