"""The image of the two matrices (Karoline, 2026-10-04: "as the first result, generate and show a downloadable
image of the 2 matrices")."""
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
    circles = list(root.iter(f"{SVG}circle"))
    # presence only: B's formate and C's glucose, one open circle each, plus the one in the key
    assert len([c for c in circles if c.get("fill") == "none"]) == 3
    # every value rests on two or more replicates: no single-replicate dots
    assert len([c for c in circles if c.get("fill") != "none"]) == 0


def test_booleans_print_no_numbers(client):
    texts = _texts(figure.matrices_svg(run(client, booleans=True)))
    assert "Consumed" in texts and not any(t.replace(".", "").isdigit() for t in texts)


def test_the_page_shows_it_first_and_offers_it_as_a_download(client):
    from foodnet import gui
    page = gui.render_result("T", run(client))
    assert page.index("<h2 id=\"matrices\">Consumed and produced") < page.index("<h2 id=\"taxa\">Taxa</h2>")
    assert "/matrices.svg?token=T" in page and "Download the image (.svg)" in page


def test_the_buttons_come_first_then_an_index_of_the_sections(client):
    # Karoline, 2026-10-07: "buttons related to results should appear above the matrix images and there should
    # be an index below the buttons to allow users to jump to different result sections"
    import re

    from foodnet import gui
    page = gui.render_result("T", run(client, report_rates=True))
    result = page[page.index('id="result"'):]
    buttons, nav, figure = (result.index("Download network"), result.index('<nav class="index"'),
                            result.index('id="matrices"'))
    assert buttons < nav < figure and result.index("Get CRM parameters") < nav
    links = re.findall(r'<a href="#([a-z]+)">', result[nav:result.index("</nav>", nav)])
    assert links[:4] == ["matrices", "taxa", "arcs", "sources"]
    assert all(f'id="{anchor}"' in result for anchor in links)


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


def test_every_cell_names_its_source_studies_on_mouseover(client):
    # Karoline, 2026-10-06: "make fields in the output matrices interactive, so a mouseover will show the
    # source studies". The text is the cell's aria-label and its hover box, never an SVG title as well (2026-10-07:
    # "hover above the matrix images now produces 2 boxes instead of one").
    root = ET.fromstring(figure.matrices_svg(run(client)))
    assert not list(root.iter(f"{SVG}title"))
    titles = [g.get("aria-label") for g in root.iter(f"{SVG}g") if g.get("class") == "fn-cell"]
    acetate = next(t for t in titles if t.startswith("Alpha alpha A1, acetate, produced"))
    assert "4 \u00b1 0.1 mM, 3 replicate(s)" in acetate          # both replicates +4, so sd 0
    assert "SMGDB00000001 (Synthetic study one)" in acetate and "EMGDB000000001" in acetate
    formate = next(t for t in titles if t.startswith("Beta beta B1, formate, produced"))
    assert "not assayed" in formate and "Seen in mMCB (+5 mM, 2 replicate(s)): SMGDB00000002" in formate
    # one per cell: 3 taxa x 4 metabolites x 2 matrices
    assert len(titles) == 24


def test_a_value_from_another_medium_has_its_own_background(client):
    root = ET.fromstring(figure.matrices_svg(run(client, presence_entries="value")))
    fills = [r.get("fill") for r in root.iter(f"{SVG}rect")]
    assert fills.count(figure.OTHER_MEDIUM) == 3       # B's formate, C's glucose, and the key
    texts = _texts(figure.matrices_svg(run(client, presence_entries="value")))
    assert "5.0" in texts and "value from another medium, not comparable with the gray scale" in texts


def test_the_mouseover_text_is_drawn_in_the_image_not_left_to_the_browser(client):
    # Karoline, 2026-10-06: "the hover didn't land": a browser need not show an SVG title, so each cell also
    # has a box with the same text, hidden until the cell is hovered, by the image's own style sheet
    svg = figure.matrices_svg(run(client))
    root = ET.fromstring(svg)
    tips = [g for g in root.iter(f"{SVG}g") if g.get("class") == "fn-tip"]
    assert len(tips) == 24 and all(g.get("visibility") == "hidden" for g in tips)
    assert "#fn-c1:hover ~ #fn-t1" in svg and "visibility: visible" in svg
    first = " ".join("".join(t.itertext()) for t in tips[0].iter(f"{SVG}text"))
    assert first.startswith("Alpha alpha A1, acetate, consumed")
