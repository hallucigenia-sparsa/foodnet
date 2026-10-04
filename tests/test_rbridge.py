"""Sending CRM parameters to a listening R session (foodnet.rbridge).

The wire and the answer are checked here, with a stand-in listener that answers the way the R package does;
what the R package does with the payload is checked in r/tests/testthat.
"""
import http.server
import json
import threading

import pytest
from conftest import run

from foodnet import matrix, rbridge


class _Listener(http.server.BaseHTTPRequestHandler):
    """What `foodnet::foodnet_listen` does, in Python: take the POST and answer what it holds."""

    received: list = []
    answer = {"received": True, "taxa": 3, "resources": 4, "growth_rates": 1, "without_a_rate": 2}

    def log_message(self, *_args):
        pass

    def do_POST(self):            # noqa: N802 - the name http.server requires
        body = self.rfile.read(int(self.headers.get("Content-Length") or 0))
        type(self).received.append((self.path, json.loads(body.decode("utf-8"))))
        payload = json.dumps(self.answer).encode("utf-8")
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(payload)))
        self.end_headers()
        self.wfile.write(payload)


@pytest.fixture
def listener():
    """A stand-in R session on a free port; yields (port, what it received)."""
    _Listener.received = []
    server = http.server.HTTPServer(("127.0.0.1", 0), _Listener)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    yield server.server_address[1], _Listener.received
    server.shutdown()


def test_the_payload_reaches_r_and_r_says_what_it_holds(listener, client):
    port, received = listener
    payload = matrix.crm_payload(run(client, report_rates=True))
    answer = rbridge.send(payload, port=port)
    assert answer["received"] and answer["taxa"] == 3
    path, body = received[-1]
    assert path == "/foodnet/crm" and body["format"] == "foodnet.crm/v0"


def test_nothing_listening_says_what_to_run_in_r():
    with pytest.raises(rbridge.RError, match="foodnet_listen"):
        rbridge.send({"format": "x"}, port=9, timeout=2)


def test_only_this_machine_is_ever_addressed():
    with pytest.raises(rbridge.RError, match="only talks to an R session on this machine"):
        rbridge._local("http://example.org:8794/foodnet/crm")


def test_the_default_port_is_not_grownets():
    assert rbridge.DEFAULT_PORT == 8794
