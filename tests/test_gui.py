"""The local page: the form, the settings it reads back, and every route through a running server."""
import io
import json
import threading
import urllib.error
import urllib.parse
import urllib.request
import zipfile

import pytest
from conftest import FakeClient

from foodnet import gui
from foodnet.search import DEFAULTS


def test_the_form_has_two_boxes_and_the_phase_choice_with_exponential_first():
    page = gui.render_form("T")
    assert 'name="species"' in page and 'name="conditions"' in page
    # Karoline, 2026-10-04: "a multiple-choice radio button below that says: 'Exponential phase' (the
    # default), 'Stationary phase', 'Both'"
    radios = [part.split('"')[0] for part in page.split('name="phase" value="')[1:]]
    assert radios == ["exponential", "stationary", "both"]
    assert 'value="exponential" checked' in page
    assert page.index('name="conditions"') < page.index('name="phase"') < page.index("<summary>Advanced settings")


def test_the_advanced_settings_hold_the_window_the_limit_and_the_media_options():
    page = gui.render_form("T")
    for name in ("window_start", "window_end", "detection_limit", "ignore_media", "booleans", "merge_arcs",
                 "merge_genera", "fraction"):
        assert f'name="{name}"' in page, name
    assert 'value="0.2"' in page


def test_settings_are_read_back_and_bad_values_fall_back_to_the_defaults():
    form = {"phase": ["both"], "window_start": ["0"], "window_end": ["48"], "fraction": ["80"],
            "detection_limit": ["0.5"], "booleans": ["1"], "min_studies": ["x"]}
    s = gui.parse_settings(form)
    assert s["phase"] == "both" and (s["window_start"], s["window_end"]) == (0.0, 48.0)
    assert s["fraction"] == pytest.approx(0.8) and s["detection_limit"] == 0.5 and s["booleans"]
    assert s["min_studies"] == DEFAULTS["min_studies"]
    # a window that ends before it starts is no window
    assert gui.parse_settings({"window_start": ["10"], "window_end": ["5"]})["window_start"] is None


def test_crm_mode_switches_growth_rates_on_and_off():
    on = gui.crm_mode(dict(DEFAULTS), True)
    assert gui.crm_mode_on(on) and not gui.crm_mode_on(gui.crm_mode(on, False))


@pytest.fixture
def server():
    handler = type("H", (gui._Handler,), {"token": "tok", "client_factory": staticmethod(FakeClient),
                                          "state": {}, "wait": 5.0})
    httpd = gui._Server(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    thread.start()
    yield f"http://127.0.0.1:{httpd.server_address[1]}"
    httpd.shutdown()


def _get(url):
    with urllib.request.urlopen(url, timeout=10) as r:
        return r.read(), r.headers


def _search(base, **extra):
    data = urllib.parse.urlencode({"species": "Alpha alpha\nBeta beta\nGamma gamma", **extra}).encode()
    with urllib.request.urlopen(f"{base}/run?token=tok", data=data, timeout=30) as r:
        return r.read().decode(), r.geturl()


def test_a_request_without_the_token_gets_the_page_that_says_where_the_link_is(server):
    with pytest.raises(urllib.error.HTTPError) as e:
        _get(f"{server}/")
    assert e.value.code == 403 and "is running on this machine" in e.value.read().decode()


def test_a_search_shows_its_arcs_and_every_download_works(server):
    page, url = _search(server)
    assert "arc(s):" in page and "presence only" in page and "Values from" in page
    job = urllib.parse.parse_qs(urllib.parse.urlparse(url).query)["job"][0]
    q = f"token=tok&job={job}"
    body, _ = _get(f"{server}/download?{q}&format=json")
    assert json.loads(body)["schema"] == "foodnet.metabolite_network/v0"
    body, _ = _get(f"{server}/download?{q}&format=graphml")
    assert body.startswith(b"<?xml")
    body, _ = _get(f"{server}/download?{q}&format=matrix")
    assert body.decode().startswith('taxon [values from Wilkins')
    body, _ = _get(f"{server}/download?{q}&format=matrices")
    assert "consumed.csv" in zipfile.ZipFile(io.BytesIO(body)).namelist()
    body, _ = _get(f"{server}/report.txt?{q}")
    assert b"foodnet" in body


def test_crm_parameters_need_crm_mode_and_then_download(server):
    page, url = _search(server)
    assert "Get CRM parameters" not in page
    page, url = _search(server, report_rates="1")
    assert "Get CRM parameters" in page and "foodnet_crm(" in page
    job = urllib.parse.parse_qs(urllib.parse.urlparse(url).query)["job"][0]
    body, _ = _get(f"{server}/crm.zip?token=tok&job={job}")
    assert "growth_rates.csv" in zipfile.ZipFile(io.BytesIO(body)).namelist()
    body, _ = _get(f"{server}/crm.json?token=tok&job={job}")
    assert json.loads(body)["format"] == "foodnet.crm/v1"


def test_help_legend_and_about_open(server):
    for path in ("help", "legend", "about"):
        body, _ = _get(f"{server}/{path}?token=tok")
        assert b"food<b>net</b>" in body


def test_the_search_button_says_what_it_gives():
    # Karoline, 2026-10-04: 'Change button label "Find metabolites" to "Get taxon-metabolite network"'
    page = gui.render_form("T")
    assert '<button class="primary" type="submit">Get taxon-metabolite network</button>' in page
    assert "Find metabolites" not in page


def test_metabolite_nodes_carry_their_chebi_id(client):
    # Karoline, 2026-10-04: "put CheBI identifiers as metabolite node attributes"
    from conftest import run
    net = run(client)["network"]
    assert {n.name: n.chebi_id for n in net.metabolites()} == {
        "acetate": "30089", "butyrate": "17968", "formate": "15740", "glucose": "17234"}


def test_a_window_that_cannot_be_used_is_said_not_dropped():
    # found in the pre-release audit: a window ending before its start was silently ignored
    assert gui.window_problem({"window_start": ["10"], "window_end": ["5"]}).startswith(
        "The time window ends at 5 h, which is not after its start at 10 h")
    assert "needs numbers" in gui.window_problem({"window_start": ["a"], "window_end": ["5"]})
    assert "needs a start" in gui.window_problem({"window_end": ["5"]})
    assert gui.window_problem({"second_window_start": ["0"], "second_window_end": [""]}) == ""
    assert gui.window_problem({"window_start": ["0"], "window_end": ["48"]}) == ""


def test_a_search_that_is_gone_is_a_404_never_another_searchs_data(server):
    _search(server)                                  # a result exists, so a fallback would have served it
    with pytest.raises(urllib.error.HTTPError) as e:
        _get(f"{server}/crm.json?token=tok&job=deadbeef")
    assert e.value.code == 404 and "no longer here" in e.value.read().decode()
    with pytest.raises(urllib.error.HTTPError) as e:
        _get(f"{server}/download?token=tok&job=deadbeef&format=json")
    assert e.value.code == 404


def test_only_web_addresses_become_study_links():
    assert gui._web("https://doi.org/10.1/x") and gui._web(" http://example.org ")
    assert not gui._web("javascript:alert(1)") and not gui._web("") and not gui._web("data:text/html,x")
