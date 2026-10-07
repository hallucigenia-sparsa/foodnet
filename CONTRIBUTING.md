# Contributing to foodnet

foodnet is a small open project of the KU Leuven Laboratory of Molecular Bacteriology, the sister tool
of grownet. Issues and pull requests are welcome. Discipline is deliberate here: the checks below run
the same way for everyone, in CI and before every commit, so the repository stays honest and portable.

## Working with coding agents

Features, tasks, and the roles agents take (mayor, worker, verifier) are described in
[AGENTS.md](AGENTS.md), with shared agent memory in [docs/agents/NOTES.md](docs/agents/NOTES.md). Humans
propose work by opening a Feature issue.

## Development setup

```
git clone https://github.com/hallucigenia-sparsa/foodnet
cd foodnet
python -m venv .venv
. .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -e ".[dev]"
pre-commit install            # optional: run the checks on every commit (pip install pre-commit first)
```

The page then opens with `python -m foodnet gui`. How a version is released (PyPI and the Windows
program) is in [RELEASING.md](RELEASING.md).

The R companion package in [`r/`](r) is built and checked with `make r-check`, where R is installed, and
CI runs the same build and `R CMD check` on every push (the `r-package` job). To
try it against a working copy, install it from the clone rather than from GitHub, since what is on a
branch is not what `install_github` reads:

```r
remotes::install_local("<this clone>/r")       # or, in a terminal: R CMD INSTALL r
```

The texts that ship with the tool, the help page and both READMEs, describe the released state, so they
give the `install_github(...)` line pinned to a release tag (`ref = "vX.Y.Z"`, set in the release commit)
or plain between releases, and never a branch of ours; `tests/test_release.py` fails if work in progress
creeps into them.

## Before you open a pull request

Run the same checks CI runs, and make sure all pass. `make check` runs all three:

```
ruff check .              # lint
python checks/gate.py     # the guardrail gate
pytest -q                 # the tests
```

If you enabled `pre-commit install`, these run automatically on every commit.

## The guardrail gate

`checks/gate.py` is self-contained (no dependency on any private path) and blocks a class of mistake per
check: committed secrets, raw or pulled experimental data (only the synthetic fixtures under
`tests/fixtures/` belong in the repository, never data pulled from mGrowthDB or shared by a collaborator),
absolute local-machine paths, imports that are not the standard library or foodnet, the house style,
and a schema contract. See [docs/DATA_GOVERNANCE.md](docs/DATA_GOVERNANCE.md) for the data rules.

House style, checked on documentation: ASCII punctuation (use a comma, colon, parentheses, or the word
"to" for a range rather than an em dash or en dash), US spelling, and direct prose with no hedging
caveats. Keep it, and the check stays quiet.

## Changing the method

How phases, values, media, pooling and presence are derived is a scientific choice. The decisions, in
the words of the people who made them, and the open questions are in
[docs/METHOD_NOTES.md](docs/METHOD_NOTES.md). Two rules hold: record what the data does not support in the
skip list rather than inventing a value, and open an issue to discuss a change before you implement it.

## The network format is a contract

`src/foodnet/model.py` and the JSON Schema at `schema/metabolite_network.schema.json` are the contract
downstream tools read. Keep it backward compatible where you can, keep every edge carrying at least one
supporting study (the edge-level attribution), and regenerate the schema file from the code if you change
it (`python -m foodnet schema --out schema/metabolite_network.schema.json`). The gate fails if the
file and the code drift, or if an emitted network does not validate against its own schema.

## Style and scope

Keep pull requests small and focused, match the surrounding code, and note user-facing changes in
[CHANGELOG.md](CHANGELOG.md).
