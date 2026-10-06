"""foodnet command line: taxon and metabolite networks from mGrowthDB batch monocultures.

  foodnet derive --taxa "Escherichia coli LF82" "Bacteroides fragilis" --out network.json
  foodnet derive --taxa Roseburia --phase both --format matrices --out roseburia.zip
  foodnet validate network.json                 # check a network against the schema
  foodnet schema --out metabolite_network.schema.json
  foodnet style --out foodnet_style.xml         # the Cytoscape style, for a downloaded GraphML
  foodnet gui                                   # the local page
"""
from __future__ import annotations

import argparse
import json
import sys

from . import rates
from .search import DEFAULTS, PHASE_LABELS

# the exit status of a derive whose result is incomplete (records that could not be read): its files are
# written and say so, but a script must not take them for a complete result
INCOMPLETE = 3

DERIVE_EXAMPLES = """examples:
  the page's Example, with its report, sent to Cytoscape:
    foodnet derive --taxa "Escherichia coli LF82" "Bacteroides fragilis" "Roseburia intestinalis" \\
        --out example.json --report example_report.txt --to-cytoscape

  both growth phases, as the consumed and produced matrices:
    foodnet derive --taxa Roseburia --phase both --format matrices --out roseburia.zip

  one time window instead of the phases, values from every medium, as booleans:
    foodnet derive --taxa Blautia Bacteroides --window 0 48 --ignore-media --booleans --format matrix \\
        --out blautia_bacteroides.csv

  values from Wilkins-Chalgren only, with growth rates, as CRM parameters:
    foodnet derive --taxa Blautia Roseburia --conditions "Wilkins-Chalgren" --crm-mode --crm wc_crm.zip

  the same parameters sent into a waiting R session (in R: library(foodnet); foodnet_listen()):
    foodnet derive --taxa Blautia Roseburia --crm-mode --to-r

"""


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="foodnet", description=__doc__.split("\n\n")[0])
    sub = p.add_subparsers(dest="command")

    d = sub.add_parser("derive", help="derive a network from mGrowthDB (live)", epilog=DERIVE_EXAMPLES,
                       formatter_class=argparse.RawDescriptionHelpFormatter)
    d.add_argument("--taxa", nargs="+", default=[], help="species, strains, genera or NCBI taxon ids")
    d.add_argument("--all", dest="all_studies", action="store_true",
                   help="every batch monoculture with metabolites in mGrowthDB (the page's All)")
    d.add_argument("--conditions", default="", help="media, experiments or studies to limit the search to "
                   "(comma separated); empty: all data, values from the medium holding data for the most taxa")
    d.add_argument("--phase", choices=list(PHASE_LABELS), default=DEFAULTS["phase"])
    d.add_argument("--window", nargs=2, type=float, metavar=("START", "END"),
                   help="a time window in hours that replaces the phases")
    d.add_argument("--second-window", nargs="+", metavar="NAME",
                   help="metabolites measured over a second time window (see --second-window-from/--to)")
    d.add_argument("--second-window-from", type=float, default=0.0, help="its start in hours (default 0)")
    d.add_argument("--second-window-to", type=float, default=None,
                   help="its end in hours (default: each culture's last sample)")
    d.add_argument("--fraction", type=float, default=DEFAULTS["fraction"],
                   help="exponential growth ends at this share of the maximal abundance (default 0.9)")
    d.add_argument("--no-growth-factor", type=float, default=DEFAULTS["no_growth_factor"])
    d.add_argument("--detection-limit", type=float, default=DEFAULTS["detection_limit"], help="mM (default 0.2)")
    d.add_argument("--mean-only", action="store_true",
                   help="judge a change by its mean alone, not also by the replicates' spread")
    d.add_argument("--compound-limits", default=DEFAULTS["compound_limits"],
                   help="detection limits of their own, as name=mM, comma separated: thiamine=0.01")
    d.add_argument("--ignore-media", action="store_true", help="values from every medium, pooled")
    d.add_argument("--outside-evidence", action="store_true",
                   help="with --conditions, also take presence-only evidence from everything outside it")
    d.add_argument("--booleans", action="store_true", help="report 1, 0 or NA instead of amounts")
    d.add_argument("--presence-entries", choices=["na", "true", "value"], default=DEFAULTS["presence_entries"],
                   help="how a cell seen only in another medium is written in the matrices (default na)")
    d.add_argument("--report-rates", action="store_true", help="collect growth rates")
    d.add_argument("--crm-mode", action="store_true", help="what a consumer-resource model needs (growth rates on)")
    d.add_argument("--rate-method", choices=list(rates.METHODS), default=DEFAULTS["rate_method"])
    d.add_argument("--rate-window", type=int, default=DEFAULTS["rate_window"])
    d.add_argument("--merge-arcs", action="store_true", help="one arc across studies")
    d.add_argument("--min-studies", type=int, default=DEFAULTS["min_studies"])
    d.add_argument("--merge-genera", action="store_true", help="one node per genus")
    d.add_argument("--exclude-studies", default=DEFAULTS["exclude_studies"])
    d.add_argument("--exclude-experiments", default=DEFAULTS["exclude_experiments"],
                   help="experiment ids never used, comma separated")
    d.add_argument("--no-strict-media", action="store_true",
                   help="tell media apart by their names only, not by altered composition or atmosphere")
    d.add_argument("--exclude-metabolites", default=DEFAULTS["exclude_metabolites"])
    d.add_argument("--include-non-batch", action="store_true",
                   help="also read chemostat and serial dilution monocultures")
    d.add_argument("--spike-factor", type=float, default=DEFAULTS["spike_factor"])
    d.add_argument("--correction", choices=["bh", "by"], default=DEFAULTS["correction"])
    d.add_argument("--format", choices=["json", "graphml", "matrix", "matrices"], default="json")
    d.add_argument("--out", help="where to write the network (default: standard output; matrices need a file)")
    d.add_argument("--report", help="also write the report to this file")
    d.add_argument("--figure", help="also write the consumed and produced matrices as an image (SVG)")
    d.add_argument("--rates", help="also write the growth rates (CSV); needs --report-rates or --crm-mode")
    d.add_argument("--crm", help="also write the CRM parameters (zip); needs --report-rates or --crm-mode")
    d.add_argument("--to-r", action="store_true", help="send the CRM parameters to a listening R session")
    d.add_argument("--r-port", type=int, default=None)
    d.add_argument("--allow-incomplete", action="store_true",
                   help=f"exit 0 even when mGrowthDB records could not be read (otherwise {INCOMPLETE}; the files "
                        "say incomplete either way)")
    d.add_argument("--to-cytoscape", action="store_true", help="send the network to a running Cytoscape")
    d.add_argument("--cytoscape-port", type=int, default=1234)

    v = sub.add_parser("validate", help="check a network file against the schema")
    v.add_argument("file")
    s = sub.add_parser("schema", help="write the network JSON Schema")
    s.add_argument("--out")
    st = sub.add_parser("style", help="write the Cytoscape style (XML)")
    st.add_argument("--out", default="foodnet_style.xml")
    g = sub.add_parser("gui", help="open the local page")
    g.add_argument("--port", type=int, default=0)
    g.add_argument("--no-browser", action="store_true")
    return p


def settings_from(a) -> dict:
    s = dict(DEFAULTS)
    s.update(phase=a.phase, fraction=a.fraction, no_growth_factor=a.no_growth_factor,
             detection_limit=a.detection_limit, ignore_media=a.ignore_media, booleans=a.booleans,
             judge_spread=not a.mean_only, compound_limits=a.compound_limits,
             outside_evidence=a.outside_evidence, presence_entries=a.presence_entries,
             second_window_metabolites=", ".join(a.second_window or []),
             second_window_start=a.second_window_from, second_window_end=a.second_window_to,
             report_rates=a.report_rates or a.crm_mode, rate_method=a.rate_method, rate_window=a.rate_window,
             merge_arcs=a.merge_arcs, min_studies=a.min_studies, merge_genera=a.merge_genera,
             conditions="\n".join(x.strip() for x in a.conditions.split(",") if x.strip()),
             exclude_studies=a.exclude_studies, exclude_metabolites=a.exclude_metabolites,
             exclude_experiments=a.exclude_experiments, strict_media=not a.no_strict_media,
             include_non_batch=a.include_non_batch, spike_factor=a.spike_factor, correction=a.correction)
    from .search import compound_limits_of
    try:
        compound_limits_of(s)
    except ValueError as e:
        raise SystemExit(f"--compound-limits: {e}") from None
    if a.window:
        start, end = a.window
        if end <= start:
            raise SystemExit("--window END must be after START")
        s["window_start"], s["window_end"] = start, end
    return s


def _write(path, data) -> None:
    mode = "wb" if isinstance(data, bytes) else "w"
    with open(path, mode, **({} if isinstance(data, bytes) else {"encoding": "utf-8"})) as f:
        f.write(data)


def _derive(a) -> int:
    from . import matrix
    from .export import to_graphml
    from .mgrowthdb import MGrowthDBClient, MGrowthDBError
    from .report import report_text
    from .search import run_query

    if not a.taxa and not a.all_studies:
        print("derive needs --taxa (or --all); see foodnet derive --help", file=sys.stderr)
        return 2
    s = settings_from(a)
    if (a.rates or a.crm or a.to_r) and not s["report_rates"]:
        print("--rates, --crm and --to-r need the growth rates: add --crm-mode (or --report-rates)", file=sys.stderr)
        return 2
    if a.format == "matrices" and not a.out:
        print("--format matrices writes a zip, so it needs --out", file=sys.stderr)
        return 2

    def progress(done, total, message):
        print(f"  {message}", file=sys.stderr)

    try:
        result = run_query(MGrowthDBClient(), a.taxa, s, progress=progress, all_studies=a.all_studies)
    except MGrowthDBError as e:
        print(f"mGrowthDB is not reachable: {e}", file=sys.stderr)
        return 1
    for line in result["errors"] + result["warnings"]:
        print(f"note: {line}", file=sys.stderr)
    for entry in result["unresolved"]:
        print(f"not used: {entry}: {result['reasons'].get(entry, '')}", file=sys.stderr)
    net = result["network"]
    problems = net.validate()
    if problems:
        raise SystemExit("network invalid:\n  " + "\n  ".join(problems))
    writers = {"json": lambda: net.to_json(), "graphml": lambda: to_graphml(net),
               "matrix": lambda: matrix.signed_csv(net, result), "matrices": lambda: matrix.pair_package(result)}
    out = writers[a.format]()
    if a.out:
        _write(a.out, out)
        print(f"wrote {a.out}: {len(net.taxa())} taxa, {len(net.metabolites())} metabolites, {len(net.edges)} arcs",
              file=sys.stderr)
    else:
        sys.stdout.write(out)
    if a.report:
        _write(a.report, report_text(result))
    if a.figure:
        from .figure import matrices_svg
        _write(a.figure, matrices_svg(result))
    if a.rates:
        _write(a.rates, matrix.rates_csv(result))
    if a.crm:
        _write(a.crm, matrix.crm_package(result))
    if a.to_r:
        from . import rbridge
        try:
            answer = rbridge.send(matrix.crm_payload(result), **({"port": a.r_port} if a.r_port else {}))
        except rbridge.RError as e:
            print(f"Send to R: {e}", file=sys.stderr)
            return 1
        print(f"sent to R: {answer.get('taxa', 0)} taxa, {answer.get('resources', 0)} resources", file=sys.stderr)
    if a.to_cytoscape:
        from .cytoscape import CytoscapeError, send
        try:
            sent = send(net, port=a.cytoscape_port)
        except CytoscapeError as e:
            print(f"Send to Cytoscape: {e}", file=sys.stderr)
            return 1
        print(f"sent to Cytoscape: network {sent['suid']}", file=sys.stderr)
    if result["errors"] and not a.allow_incomplete:
        # the files are written, and say so in their meta; a script must not take them for complete
        print("foodnet: the result is incomplete (mGrowthDB records could not be read; see the notes above). "
              f"Exit status {INCOMPLETE}; --allow-incomplete accepts it.", file=sys.stderr)
        return INCOMPLETE
    return 0


def main(argv=None) -> int:
    a = build_parser().parse_args(argv)
    if a.command == "derive":
        try:
            return _derive(a)
        except OSError as e:
            # a path that cannot be written is one line, not a traceback (grownet, Karoline 2026-09-28)
            print(f"foodnet: {e.filename}: {e.strerror}", file=sys.stderr)
            return 2
    if a.command == "validate":
        from .schema import validate_file
        try:
            problems = validate_file(a.file)
        except OSError as e:
            print(f"foodnet: {a.file}: {e.strerror}", file=sys.stderr)
            return 2
        except ValueError:
            print(f"foodnet: {a.file} is not valid JSON", file=sys.stderr)
            return 2
        if problems:
            print("\n".join(problems))
            return 1
        with open(a.file, encoding="utf-8") as f:
            print(f"{a.file}: valid ({json.load(f).get('schema')})")
        return 0
    if a.command == "schema":
        from .schema import schema_json
        if a.out:
            _write(a.out, schema_json())
        else:
            sys.stdout.write(schema_json())
        return 0
    if a.command == "style":
        from .cytoscape import style_xml
        _write(a.out, style_xml())
        print(f"wrote {a.out}")
        return 0
    if a.command == "gui":
        from .gui import serve
        serve(port=a.port, open_browser=not a.no_browser)
        return 0
    build_parser().print_help()
    return 0


if __name__ == "__main__":
    sys.exit(main())
