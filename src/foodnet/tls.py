"""The TLS context of foodnet's connections to mGrowthDB and ChEBI (Karoline, 2026-10-22: "go with option 1, like
grownet").

On macOS the certificates are checked against the system's trust store, through truststore, foodnet's one runtime
dependency (macOS only). Pythons that uv or the python.org installer put on a Mac check them against
/etc/ssl/cert.pem, which lacks the root mGrowthDB's certificate chains to (HARICA TLS RSA Root CA 2021), so
`uv tool install foodnet`, the README's install, failed at the first search; tests passed because the
development Python (Anaconda's) and Homebrew's ship their own bundles. Elsewhere, or where truststore is missing,
fails to load or is older than 0.10.5, it is Python's default context.

Each thread has a context of its own: truststore before 0.10.5 could leave a context shared by concurrent
handshakes with certificate checking switched off for good (a review measured it on 0.10.0 and 0.10.4; its
changelog: "Fixed a thread-safety issue for concurrent handshakes on macOS and Windows"), and even 0.10.5
holds one lock over every handshake of a shared context.
"""
from __future__ import annotations

import ssl
import sys
import threading

TRUSTSTORE_MIN = (0, 10, 5)
_local = threading.local()


def _version(text: str) -> tuple:
    parts = []
    for part in text.split(".")[:3]:
        digits = "".join(ch for ch in part if ch.isdigit())
        parts.append(int(digits) if digits else 0)
    return tuple(parts)


def _make() -> ssl.SSLContext:
    if sys.platform == "darwin":
        try:
            import importlib.metadata

            import truststore
            if _version(importlib.metadata.version("truststore")) >= TRUSTSTORE_MIN:
                return truststore.SSLContext(ssl.PROTOCOL_TLS_CLIENT)
        except Exception:  # noqa: BLE001 - missing, too old, or failing to load: the default context still verifies
            pass
    return ssl.create_default_context()


def context() -> ssl.SSLContext:
    """This thread's context."""
    made = getattr(_local, "context", None)
    if made is None:
        made = _local.context = _make()
    return made


def system_store() -> bool:
    """Whether certificates are checked against the operating system's store (macOS with truststore)."""
    return type(context()).__module__.split(".")[0] == "truststore"
