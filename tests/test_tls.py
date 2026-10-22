"""foodnet's TLS context: the macOS trust store through truststore, else Python's default (Karoline, 2026-10-22)."""
import ssl
import sys
import types

from foodnet import tls


def _fresh(monkeypatch, platform, truststore, version="0.10.5"):
    import importlib.metadata
    monkeypatch.setattr(tls, "_local", __import__("threading").local())
    monkeypatch.setattr(sys, "platform", platform)
    monkeypatch.setitem(sys.modules, "truststore", truststore)
    monkeypatch.setattr(importlib.metadata, "version", lambda name: version)


def _fake_truststore(made):
    class Context(ssl.SSLContext):
        pass
    Context.__module__ = "truststore._api"
    return types.SimpleNamespace(SSLContext=lambda protocol: made.append(protocol) or Context(protocol))


def test_macos_checks_certificates_against_the_system_store(monkeypatch):
    made = []
    _fresh(monkeypatch, "darwin", _fake_truststore(made))
    assert tls.system_store() and made == [ssl.PROTOCOL_TLS_CLIENT]


def test_without_truststore_too_old_or_off_macos_it_is_pythons_default(monkeypatch):
    # a review: truststore before 0.10.5 could switch certificate checking off under concurrent handshakes
    for platform, module, version in (("darwin", None, "0.10.5"), ("linux", _fake_truststore([]), "0.10.5"),
                                      ("darwin", _fake_truststore([]), "0.10.4"),
                                      ("darwin", types.SimpleNamespace(), "0.10.5")):   # fails to load
        _fresh(monkeypatch, platform, module, version)     # None in sys.modules makes the import fail
        context = tls.context()
        assert type(context) is ssl.SSLContext and context.verify_mode == ssl.CERT_REQUIRED
        assert context.check_hostname and not tls.system_store()


def test_each_thread_has_a_context_of_its_own(monkeypatch):
    import threading
    _fresh(monkeypatch, "darwin", _fake_truststore([]))
    seen = []
    threads = [threading.Thread(target=lambda: seen.append(tls.context())) for _ in range(4)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    assert len({id(c) for c in seen}) == 4 and tls.context() is tls.context()


def test_mgrowthdb_connections_use_that_context():
    from foodnet.mgrowthdb import MGrowthDBClient
    client = MGrowthDBClient()
    conn = client._connection()
    try:
        assert conn._context is tls.context()
    finally:
        conn.close()


def test_pyproject_declares_what_the_gate_allows():
    import pathlib
    import re

    sys.path.insert(0, str(pathlib.Path(__file__).parents[1] / "checks"))
    import gate
    text = (pathlib.Path(__file__).parents[1] / "pyproject.toml").read_text()
    declared = re.search(r"^dependencies = \[(.*)\]$", text, re.M).group(1)
    names = {re.match(r"[A-Za-z0-9_.-]+", d.strip().strip("\"'")).group(0)
             for d in re.findall(r"[\"']([^\"']+)[\"']", declared)}
    assert names == set(gate.DECLARED)
    assert "sys_platform == 'darwin'" in declared
