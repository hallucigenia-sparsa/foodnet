"""foodnet's TLS context: the macOS trust store through truststore, else Python's default (Karoline, 2026-10-22)."""
import ssl
import sys
import types

from foodnet import tls


def _fresh(monkeypatch, platform, truststore):
    tls.context.cache_clear()
    monkeypatch.setattr(sys, "platform", platform)
    monkeypatch.setitem(sys.modules, "truststore", truststore)


def test_macos_checks_certificates_against_the_system_store(monkeypatch):
    made = []

    class Context(ssl.SSLContext):
        pass
    Context.__module__ = "truststore._api"
    fake = types.SimpleNamespace(SSLContext=lambda protocol: made.append(protocol) or Context(protocol))
    _fresh(monkeypatch, "darwin", fake)
    try:
        assert tls.system_store() and made == [ssl.PROTOCOL_TLS_CLIENT]
    finally:
        tls.context.cache_clear()


def test_without_truststore_or_off_macos_it_is_pythons_default(monkeypatch):
    for platform, module in (("darwin", None), ("linux", types.SimpleNamespace(SSLContext=None))):
        _fresh(monkeypatch, platform, module)        # None in sys.modules makes the import fail
        try:
            context = tls.context()
            assert type(context) is ssl.SSLContext and context.verify_mode == ssl.CERT_REQUIRED
            assert not tls.system_store()
        finally:
            tls.context.cache_clear()


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
    names = {re.match(r"[A-Za-z0-9_.-]+", d.strip().strip("\"'")).group(0) for d in declared.split(",") if d.strip()}
    assert names == gate.DECLARED
    assert "sys_platform == 'darwin'" in declared
