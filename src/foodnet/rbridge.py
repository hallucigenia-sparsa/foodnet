"""Send consumer-resource model parameters into a running R session.

grownet's R bridge, carrying CRM parameters instead of gLV ones (Karoline, 2026-10-04: "this time the
purpose is to parameterize a consumer-resource model (CRM), also with miaSim ... as an example simulator").
The R session opens a small port on this machine, foodnet POSTs the payload (`foodnet.matrix.crm_payload`)
to it, and nothing is hosted and nothing leaves the machine. The companion package is in `r/` of this
repository; in R, `foodnet::foodnet_listen()` is what answers here.

When nothing answers, the error says the other way round: the R package can fetch the same payload from
the page's own URL (`foodnet::foodnet_crm(url)`), so a user who cannot open a port is not stuck.
"""
from __future__ import annotations

import http.client
import json
import urllib.error
import urllib.request

from .brand import REPOSITORY_SLUG

# The port the R package listens on by default. Not 1234 (Cytoscape), not grownet's 8793 (so both R
# packages can listen at once) and not the page's own port; both sides take the number as a setting.
DEFAULT_PORT = 8794
PATH = "/foodnet/crm"
# where the companion package lives and how it is installed, in one place, so the page, the help and the
# error below say the same line
# What a reader of a release should run. The help ships with the tool, so it says what holds for the
# released version and never for the branch it was written on (Karoline, 2026-10-04: "The help should
# refer to the stage the tool is in when released"). Installing from an unmerged branch is a development
# step; it is in docs/agents/NOTES.md and in CONTRIBUTING.md, not here.
INSTALL_R = f'remotes::install_github("{REPOSITORY_SLUG}", subdir = "r")'
# Why an install can fail on a public repository: a GitHub token stored on the machine, for another
# account or scope, makes GitHub answer 404 instead of serving it anonymously. The way around it is a
# local install, which needs no GitHub access; clearing the token is not suggested, since that changes
# the session for everything else in it (Karoline, 2026-10-04: "I don't think we should recommend it for
# users, as it alters their system settings in ways that can affect them negatively").
INSTALL_TROUBLE = ('If that fails with "HTTP error 404", a GitHub token stored on this machine is being '
                   'used and cannot see the repository. Installing from a clone needs no GitHub access: '
                   'remotes::install_local("<the repository>/r"), or R CMD INSTALL r in a terminal.')



class RError(RuntimeError):
    """R could not be reached or refused the parameters, with what to do about it."""


# What can go wrong on the wire: urllib's own errors, http.client's for a listener that answers with
# something that is not HTTP, and the socket errors outside urllib's wrapping (as in `cytoscape`).
WIRE_ERRORS = (urllib.error.URLError, http.client.HTTPException, TimeoutError, OSError)


def base_url(port: int = DEFAULT_PORT) -> str:
    return f"http://127.0.0.1:{int(port)}"


def _local(url: str) -> str:
    """Every request goes to this machine, checked where the request is made."""
    if not url.startswith("http://127.0.0.1:"):
        raise RError(f"refusing to send to {url!r}: foodnet only talks to an R session on this machine")
    return url


def send(payload: dict, port: int = DEFAULT_PORT, timeout: float = 30.0) -> dict:
    """POST the CRM payload to a listening R session and return what it answered.

    The answer is R's own: {"received": true, "taxa": n, "growth_rates": n, ...} from
    `foodnet::foodnet_listen`. Raises `RError` with what to do when nothing is listening or the
    listener is not the R package.
    """
    url = _local(f"{base_url(port)}{PATH}")
    body = json.dumps(payload).encode("utf-8")
    request = urllib.request.Request(url, data=body, method="POST",
                                     headers={"Content-Type": "application/json"})
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:   # noqa: S310 - localhost only
            text = response.read().decode("utf-8")
    except urllib.error.HTTPError as e:
        raise RError(f"the R session refused the parameters ({e.code} {e.reason}). Check that the foodnet "
                     "package there is up to date.") from None
    except http.client.HTTPException as e:
        raise RError(f"something is listening on port {port}, but it did not answer the way the foodnet R "
                     f"package does ({type(e).__name__}). Check which port foodnet_listen() printed and "
                     "that no other program holds it.") from None
    except WIRE_ERRORS as e:
        raise RError(unreachable(port, getattr(e, "reason", e))) from None
    try:
        answer = json.loads(text) if text.strip() else {}
    except ValueError:
        raise RError(f"the listener on port {port} answered something that is not JSON, so it is not the "
                     "foodnet R package") from None
    if not isinstance(answer, dict) or not answer.get("received"):
        raise RError(f"the listener on port {port} did not confirm the parameters: {text[:200]!r}")
    return answer


def unreachable(port: int = DEFAULT_PORT, reason: object = "") -> str:
    """What to do when nothing answers, worded for the page and for the command line."""
    detail = f" ({reason})" if reason else ""
    return (f"no R session is listening on port {port}{detail}. In R: "
            f"install.packages(\"remotes\"); {INSTALL_R}; library(foodnet); crm <- foodnet_listen(). "
            "Or fetch the same parameters from R without a listener, with foodnet_crm(url), where url is "
            "the CRM address shown under this control.")
