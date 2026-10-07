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
import os
import platform
import re
import urllib.error
import urllib.request

from . import __version__
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
#
# A release pins the line to its own tag, so the R package installed is the one this version speaks to: the
# listener's secret came in 0.2.0, and an unpinned line gave 0.1.0 users a newer R package that refused
# their page (a review). A development version names no tag, since it has none.
INSTALL_R = (f'remotes::install_github("{REPOSITORY_SLUG}", subdir = "r", ref = "v{__version__}")'
             if re.fullmatch(r"[0-9]+\.[0-9]+\.[0-9]+", __version__) else
             f'remotes::install_github("{REPOSITORY_SLUG}", subdir = "r")')
# Why an install can fail on a public repository: a GitHub token stored on the machine, for another
# account or scope, makes GitHub answer 404 instead of serving it anonymously. The way around it is a
# local install, which needs no GitHub access; clearing the token is not suggested, since that changes
# the session for everything else in it (Karoline, 2026-10-04: "I don't think we should recommend it for
# users, as it alters their system settings in ways that can affect them negatively").
INSTALL_TROUBLE = ('If that fails with "HTTP error 404" or "cannot open URL", either the release it names is not '
                   'published yet, or a GitHub token stored on this machine is being used and cannot see the '
                   'repository. Installing from a clone needs no GitHub access: '
                   'remotes::install_local("<the repository>/r"), or R CMD INSTALL r in a terminal.')



class RError(RuntimeError):
    """R could not be reached or refused the parameters, with what to do about it."""


# What can go wrong on the wire: urllib's own errors, http.client's for a listener that answers with
# something that is not HTTP, and the socket errors outside urllib's wrapping (as in `cytoscape`).
WIRE_ERRORS = (urllib.error.URLError, http.client.HTTPException, TimeoutError, OSError)


def token_path(port: int = DEFAULT_PORT) -> str:
    """Where `foodnet_listen()` leaves the secret of its listener on `port`: R's
    `tools::R_user_dir("foodnet", "cache")`, computed the way R computes it, plus "listen-token-<port>". Base R
    cannot bind a listening port to this machine alone, so the R side accepts parameters only from a request
    that carries the secret, and only this user can read the file. R may compute another directory when its
    own startup files (~/.Renviron) set R_USER_CACHE_DIR or XDG_CACHE_HOME; FOODNET_R_TOKEN_DIR, or the
    same variables set for this program, then point here."""
    if os.environ.get("FOODNET_R_TOKEN_DIR"):
        return os.path.join(os.environ["FOODNET_R_TOKEN_DIR"], f"listen-token-{int(port)}")
    if os.environ.get("R_USER_CACHE_DIR"):
        root = os.environ["R_USER_CACHE_DIR"]
    elif os.environ.get("XDG_CACHE_HOME"):
        root = os.environ["XDG_CACHE_HOME"]
    elif os.name == "nt":
        root = os.path.join(os.environ.get("LOCALAPPDATA", ""), "R", "cache")
    elif platform.system() == "Darwin":
        root = os.path.join(os.path.expanduser("~"), "Library", "Caches", "org.R-project.R")
    else:
        root = os.path.join(os.path.expanduser("~"), ".cache")
    return os.path.join(root, "R", "foodnet", f"listen-token-{int(port)}")


def listen_token(port: int = DEFAULT_PORT, path: str | None = None) -> str:
    """The secret of the R session listening on `port`, or "" when none has written one."""
    try:
        with open(path or token_path(port), encoding="utf-8") as f:
            return f.read().strip()
    except OSError:
        return ""


def base_url(port: int = DEFAULT_PORT) -> str:
    return f"http://127.0.0.1:{int(port)}"


def _local(url: str) -> str:
    """Every request goes to this machine, checked where the request is made."""
    if not url.startswith("http://127.0.0.1:"):
        raise RError(f"refusing to send to {url!r}: foodnet only talks to an R session on this machine")
    return url


def send(payload: dict, port: int = DEFAULT_PORT, timeout: float = 60.0) -> dict:
    """POST the CRM payload to a listening R session and return what it answered.

    The answer is R's own: {"received": true, "taxa": n, "growth_rates": n, ...} from
    `foodnet::foodnet_listen`. Raises `RError` with what to do when nothing is listening or the
    listener is not the R package.
    """
    url = _local(f"{base_url(port)}{PATH}")
    body = json.dumps(payload).encode("utf-8")
    token = listen_token(port)
    if not token:
        raise RError(unreachable(port, f"no foodnet_listen() secret at {token_path(port)}; if R is listening, it "
                                       "printed where it wrote its secret: set FOODNET_R_TOKEN_DIR to that folder. "
                                       "An R package older than 0.2.0 writes none: install the one matching "
                                       f"this foodnet, {INSTALL_R}"))
    request = urllib.request.Request(url, data=body, method="POST",
                                     headers={"Content-Type": "application/json", "X-Foodnet-Token": token})
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:   # noqa: S310 - localhost only
            text = response.read().decode("utf-8")
    except urllib.error.HTTPError as e:
        if e.code == 403:
            raise RError(f"the R session on port {port} refused the parameters: the secret foodnet_listen() wrote "
                         f"was not the one sent (read from {token_path(port)}). Start foodnet_listen() again, and "
                         "check that R and this page run as the same user.") from None
        raise RError(f"the R session refused the parameters ({e.code} {e.reason}). Check that the foodnet "
                     "package there is up to date.") from None
    except http.client.HTTPException as e:
        raise RError(f"something is listening on port {port}, but it did not answer the way the foodnet R "
                     f"package does ({type(e).__name__}). Check which port foodnet_listen() printed and "
                     "that no other program holds it.") from None
    except TimeoutError:
        raise RError(f"R on port {port} did not answer within {timeout:g} s. It may be busy, or another "
                     "connection may be holding its listener; if R printed that the parameters arrived, they "
                     "did. Otherwise press Send to R again.") from None
    except WIRE_ERRORS as e:
        if isinstance(getattr(e, "reason", None), TimeoutError):
            raise RError(f"R on port {port} did not answer within {timeout:g} s. It may be busy, or another "
                         "connection may be holding its listener; if R printed that the parameters arrived, "
                         "they did. Otherwise press Send to R again.") from None
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
