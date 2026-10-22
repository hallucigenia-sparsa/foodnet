"""The TLS context of foodnet's connections to mGrowthDB and ChEBI (Karoline, 2026-10-22: "go with option 1, like
grownet").

On macOS the certificates are checked against the system's trust store, through truststore, foodnet's one runtime
dependency (macOS only). Pythons that uv or the python.org installer put on a Mac check them against
/etc/ssl/cert.pem, which lacks the root mGrowthDB's certificate chains to (HARICA TLS RSA Root CA 2021), so
`uv tool install foodnet`, the README's install, failed at the first search; tests passed because the
development Python (Anaconda's) and Homebrew's ship their own bundles. Elsewhere, or where truststore is not
installed (a clone installed without its dependencies), it is Python's default context.
"""
from __future__ import annotations

import functools
import ssl
import sys


@functools.lru_cache(maxsize=1)
def context() -> ssl.SSLContext:
    """One context for every connection (a context may be shared across threads)."""
    if sys.platform == "darwin":
        try:
            import truststore
        except ImportError:
            pass
        else:
            return truststore.SSLContext(ssl.PROTOCOL_TLS_CLIENT)
    return ssl.create_default_context()


def system_store() -> bool:
    """Whether certificates are checked against the operating system's store (macOS with truststore)."""
    return type(context()).__module__.split(".")[0] == "truststore"
