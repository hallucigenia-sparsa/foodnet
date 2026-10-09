# Releasing

Written by Claude (Karoline's agent) for the maintainers, adapted from grownet's. How a version of the tool reaches its users: on
PyPI (`uv tool install foodnet`, `pipx install foodnet`) and as a self-contained Windows program attached to
the GitHub release. The `uvx --from git+...` route in the README needs no release; it always runs `main`.

Everything is automated by `.github/workflows/release.yml`, which runs when a version tag is pushed. The CI
jobs `package` and `windows-app` build and start the same artifacts on every push, so a release only
repeats what CI has already shown to work.

## Before the first release

1. Create the GitHub repository (its address is set in `brand.REPOSITORY_SLUG` and the files named in
   docs/agents/NOTES.md) and push `main`.
2. The one-time setup below.
3. The release pull request, then the tag (below).

## One-time setup (a maintainer, in a browser)

No password or token is ever stored in the repository or in GitHub secrets: PyPI trusts the release
workflow's short-lived GitHub identity instead ("trusted publishing").

1. **PyPI account.** Sign in at <https://pypi.org> (or register), and switch on two-factor authentication,
   which PyPI requires.
2. **A pending publisher.** At <https://pypi.org/manage/account/publishing/>, under "Add a new pending
   publisher", choose GitHub and enter:
   - PyPI project name: `foodnet`
   - Owner and repository name: those of the foodnet repository
   - Workflow name: `release.yml`
   - Environment name: `pypi`

   The first release then creates the project on PyPI under that account, its first owner.
3. **A second owner.** Right after the first release, the first owner invites a second maintainer on PyPI
   (Manage project, Collaborators) with the role Owner, so the project never depends on one account.
4. **A protected environment on GitHub.** In the repository's Settings, Environments, create `pypi`. Under
   "Deployment protection rules", add the maintainers as required reviewers, so publishing waits for a
   person's approval, and under "Deployment branches and tags" allow only tags matching `v*`.

## Each release

1. **A release pull request** that:
   - sets the version in `pyproject.toml` and in `src/foodnet/__init__.py` (`__version__`);
   - turns `## [Unreleased]` in `CHANGELOG.md` into `## [x.y.z] (YYYY-MM-DD)`, and starts a new
     `## [Unreleased]` section above it if work continues;
   - sets `version` and `date-released` in `CITATION.cff` (the changelog's date);
   - sets `Version:` in `r/DESCRIPTION` to the same version (between releases it is the last release plus
     `.9000`), since the R package is installed from GitHub and its version is the only signal an installed
     copy is out of date;
   - adds `ref = "vX.Y.Z"` to the R install lines of `README.md` and `r/README.md` (the page prints the
     matching line by itself), so the R package installed is the one this version talks to; between
     releases they carry no ref, since the tag does not exist yet;
   - sets the version in README.md's links (`blob/vX.Y.Z/`, and the legend image's raw URL), since the README
     is also PyPI's page for this version, which can never be replaced;
   - adds the version, its date and a few words to `RELEASES` in `src/foodnet/help.py` (the About page's
     release history).
   `python packaging/check_release.py vX.Y.Z` must say the tag is ready (it checks every item above), and `make check` must pass:
   among other things it refuses texts that describe work in progress rather than the released state (an
   install line pointing at one of our branches, for example).
   Test the R package of the release commit before tagging with `remotes::install_local("<clone>/r")`: its
   pinned install line names a tag that exists only once step 2 is done.
2. **Merge it, then tag `main`:**

   ```bash
   git switch main && git pull
   git tag vX.Y.Z
   git push origin vX.Y.Z
   ```

3. **Watch the workflow** (Actions, "release"). It checks the tag, builds and tests the wheel on Linux,
   Windows and macOS, builds and starts the Windows program, then waits for approval of the `pypi`
   environment. Approve it; it publishes to PyPI and creates the GitHub release, with the notes taken from
   the changelog and the wheel, the source archive and `foodnet-vX.Y.Z-windows.zip` attached.
4. **Check as a user would:** `uv tool install foodnet` then `foodnet gui` on a clean machine, and the zip on
   a Windows machine. A version on PyPI cannot be replaced: a mistake is fixed with a new version.

## After the first release: signing the Windows program

The program is unsigned, so Windows warns on first start (the zip's README and the README explain it
beforehand). The route agreed on #26 and #63:

1. **The SignPath Foundation** signs open-source projects free of charge and requires an existing release
   in the form to be signed, which the first release provides. Apply at <https://signpath.org>. **Not yet
   (Karoline, 2026-09-29):** the application asks the project to show that it "is widely used or trusted"
   (media coverage, blog posts, download statistics, GitHub insights, community discussions), which a
   project released that day cannot; "If we don't have that by the next release, we should explain to
   users how to click through the warning instead". So each release's notes open with how to get past the
   warning (`packaging/check_release.py`), as the README and the zip's README.txt do. Also needed when
   applying: a public code signing policy page (the SignPath attribution line, the team by role, and a
   privacy statement saying what the program sends where), two-factor authentication for every role,
   and a product name and version in `foodnet.exe`'s file metadata. Signing
   lets Windows build up reputation across releases (an unsigned program starts from zero every time),
   should let it pass Smart App Control on Windows 11, and reduces antivirus false alarms. The publisher shown to
   a user is then the SignPath Foundation. The workflow gains a signing step before the zip is made.
2. **The Microsoft Store** (free, MSIX packaging) removes the warning entirely; worth it once the method has
   settled.
3. A commercial certificate is never bought: it would not remove the warning.
