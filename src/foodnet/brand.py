"""The tool's name, mark and page palette, in one place.

foodnet is the sister tool of grownet and shares its page layout, with its own colors so the two are never
mistaken for each other (Karoline, 2026-10-04: "can look similar in style, but with a different color
scheme").

Two signal colors, each meaning one thing everywhere (arcs, legend, Cytoscape, the page): a **produced**
arc is blue and a **consumed** arc amber. Blue against orange is the pair that stays apart under the
common color vision deficiencies (`tests/test_palette.py` simulates them), and neither is grownet's green or
orange-red. The page's primary action takes the produced blue, as grownet's takes its green.
"""
from __future__ import annotations

import colorsys
import re

NAME = "foodnet"
COMMAND = "foodnet"            # the installed command
# where the repository lives, in one place: the help, the R install line and the packaging read it here
REPOSITORY_SLUG = "hallucigenia-sparsa/foodnet"
REPOSITORY = f"https://github.com/{REPOSITORY_SLUG}"

INK, MUTED, LINE, PANEL, PAGE = "#1F1D1A", "#615B53", "#DDD8D0", "#FAF8F4", "#F2EFE9"
PRODUCED, CONSUMED = "#2160A8", "#B45309"
TAXON_NODE, METABOLITE_NODE = "#7D7468", "#E7E1D6"

_GROWNET_GENUS_COLORS = (
    "#2A78D6", "#EDA100", "#E87BA4", "#4A3AA7", "#A2C5FF", "#762E61", "#069CE4", "#C7CA85", "#B8892D",
    "#635A93", "#F1ACCC", "#A3B472", "#A74FBB", "#8AACE4", "#583A84", "#D75EB4", "#854E73", "#B25977",
    "#CC96C6", "#544EC5", "#9DDA4F", "#738242", "#8080FC", "#3C561C", "#92689C", "#4CDBE3", "#88194A",
    "#878CC9", "#4C701A", "#84C030", "#E5598E", "#BA93FB", "#959754", "#6B2094", "#724AAB", "#984260",
    "#9E658B", "#7B67CC", "#7E0F7A", "#7A5283", "#C8729C", "#E38AB5", "#B249AC", "#8E35A1", "#64B8D2",
    "#6C4302", "#A76C12", "#623B6B"
)


def _away_from_arcs(color: str) -> bool:
    """A genus color that does not read as one of the two arc colors: grownet's list, without the blues and
    ambers (hue within 25 degrees of either arc color, unless nearly gray)."""
    r, g, b = (int(color[i:i + 2], 16) / 255 for i in (1, 3, 5))
    h, _, s = colorsys.rgb_to_hls(r, g, b)
    if s < 0.25:
        return True
    hue = h * 360
    for arc in (PRODUCED, CONSUMED):
        ar, ag, ab = (int(arc[i:i + 2], 16) / 255 for i in (1, 3, 5))
        arc_hue = colorsys.rgb_to_hls(ar, ag, ab)[0] * 360
        if min(abs(hue - arc_hue), 360 - abs(hue - arc_hue)) < 25:
            return False
    return True


# Node colors by genus in Cytoscape, one per genus, taken from grownet's list (picked greedily for
# separation, most distinct first) and kept in that order, minus the colors that would read as an arc.
GENUS_COLORS = tuple(c for c in _GROWNET_GENUS_COLORS if _away_from_arcs(c))

# The mark (docs/logo.svg): a taxon (circle) produces a metabolite (square), which another taxon consumes.
# The same geometry as grownet's mark, so the two read as siblings.
LOGO = ('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 64 64" width="{size}" height="{size}" '
        'role="img" aria-label="foodnet">'
        '<g>'
        f'<line x1="16.6" y1="42.5" x2="24.9" y2="27.8" stroke="{PRODUCED}" stroke-width="5.0" stroke-linecap="round"/>'
        f'<path d="M 27.80 22.51 L 27.33 34.84 L 17.55 29.38 z" fill="{PRODUCED}" '
        f'stroke="{PRODUCED}" stroke-width="0.6" stroke-linejoin="round"/>'
        f'<line x1="36.2" y1="21.1" x2="43.7" y2="31.9" stroke="{CONSUMED}" stroke-width="5.0" stroke-linecap="round"/>'
        f'<path d="M 47.12 36.92 L 36.26 31.04 L 45.48 24.69 z" fill="{CONSUMED}" '
        f'stroke="{CONSUMED}" stroke-width="0.6" stroke-linejoin="round"/>'
        f'<circle cx="13" cy="49" r="6.4" fill="{TAXON_NODE}"/>'
        f'<rect x="25.8" y="8.8" width="12.4" height="12.4" rx="2" fill="{TAXON_NODE}"/>'
        f'<circle cx="52" cy="44" r="6.4" fill="{TAXON_NODE}"/></g>'
        '</svg>')

# the wordmark: all lowercase in the system font, with "net" in the produced blue
WORDMARK = '<span class="word">food<b>net</b></span>'

# The name in running text, marked as a name as grownet does (Karoline, 2026-10-03, for grownet: "please use
# a special style for grownet, so sentences starting with it don't look strange").
NAME_HTML = '<span class="name">food<b>net</b></span>'
# elements whose text is a command, a path, a document title or a link label: the name stays plain there.
# A link is already colored, so the name's own green inside one reads as a smudge rather than a name.
_LITERAL = ("code", "pre", "title", "script", "style", "textarea", "option", "a", "svg")
# Both patterns are bounded so that neither can be made slow by a long run of "<" or of spaces: a tag
# cannot contain another "<", and the whitespace a tag may carry before its name is a few characters at
# most (CodeQL flagged the unbounded forms as polynomial on uncontrolled data, 2026-10-04). The page
# escapes everything a user types, so no "<" of theirs reaches this, and bounding it costs nothing.
_PIECES = re.compile(r"(<[^<>]*>)")
_TAG_NAME = re.compile(r"<\s{0,8}(/?)\s{0,8}([a-zA-Z0-9]+)")
_THE_NAME = re.compile(r"(?<![\w/.-])" + NAME + r"(?![\w/.-])")


def in_prose(html: str) -> str:
    """`html` with the name styled wherever it stands as a word in a sentence.

    Only the text between tags is touched, and not inside the elements that hold commands, paths or the
    document title, so `foodnet gui` in a code block and the page's own title stay as they are.
    """
    out, literal = [], []
    for piece in _PIECES.split(html):
        tag = _TAG_NAME.match(piece)
        if tag:
            closing, name = tag.group(1), tag.group(2).lower()
            if name in _LITERAL:
                if closing:
                    if literal and literal[-1] == name:
                        literal.pop()
                elif not piece.rstrip().endswith("/>"):
                    literal.append(name)
            out.append(piece)
        elif literal:
            out.append(piece)
        else:
            out.append(_THE_NAME.sub(NAME_HTML, piece))
    return "".join(out)


def logo_svg(size: int = 64) -> str:
    """The mark at `size` pixels."""
    return LOGO.replace("{size}", str(size))


CSS = f"""
:root {{ --ink: {INK}; --muted: {MUTED}; --line: {LINE}; --panel: {PANEL}; --made: {PRODUCED};
        --took: {CONSUMED}; }}
* {{ box-sizing: border-box; }}
body {{ font: 16px/1.55 system-ui, -apple-system, "Segoe UI", sans-serif; color: var(--ink);
       background: {PAGE}; margin: 0; padding: 1.5rem 1rem; }}
.app {{ max-width: 72rem; margin: 0 auto; background: #fff; border: 1px solid var(--line);
       border-radius: 10px; overflow: hidden; }}
.app > header {{ display: flex; align-items: center; gap: .6rem; flex-wrap: wrap; padding: .8rem 1.2rem;
                border-bottom: 1px solid var(--line); background: var(--panel); }}
.brand {{ display: flex; align-items: center; gap: .55rem; margin: 0; font-size: 1.15rem;
         color: inherit; text-decoration: none; }}
.word {{ font-weight: 650; letter-spacing: -.01em; }} .word b {{ font-weight: 650; color: var(--made); }}
.name {{ font-weight: 600; }} .name b {{ font-weight: 600; color: var(--made); }}
.version {{ font-size: .8rem; color: var(--muted); border: 1px solid var(--line); border-radius: 999px;
           padding: .05rem .5rem; background: #fff; }}
.app > header nav {{ margin-left: auto; display: flex; gap: .4rem; }}
.app > main {{ padding: 1.2rem; }}
h2 {{ font-size: 1.05rem; margin: 1.8rem 0 .5rem; }}
h2.page {{ font-size: 1.3rem; margin-top: 0; }}
label.field {{ display: block; font-weight: 600; margin-bottom: .1rem; }}
.examples {{ color: var(--muted); font-size: .82rem; margin: 0 0 .4rem; }}
textarea, input, select {{ font: inherit; }}
textarea {{ width: 100%; padding: .6rem .7rem; border: 1px solid var(--line); border-radius: 8px; resize: vertical; }}
input[type=text] {{ padding: .25rem .45rem; border: 1px solid var(--line); border-radius: 6px; }}
select {{ padding: .25rem .45rem; border: 1px solid var(--line); border-radius: 6px; background: #fff; }}
textarea:focus, input:focus, select:focus, button:focus, .btn:focus {{ outline: 2px solid var(--made);
  outline-offset: 2px; }}
.hint, .muted {{ color: var(--muted); font-size: .88rem; }} .hint {{ margin: .35rem 0 0; }}
.btn, button {{ font: inherit; display: inline-block; padding: .45rem 1rem; border-radius: 8px;
               border: 1px solid var(--line); background: #fff; color: var(--ink); cursor: pointer;
               text-decoration: none; }}
.btn.primary, button.primary {{ background: var(--made); border-color: var(--made); color: #fff; font-weight: 600; }}
.btn.quiet {{ background: transparent; padding: .3rem .7rem; }}
.bar {{ display: flex; gap: .6rem; margin-top: .9rem; flex-wrap: wrap; align-items: center; }}
.beside {{ display: inline-flex; gap: .35rem; align-items: center; white-space: nowrap; }}
/* a mode is a switch: the track is green when it is on, white when the settings are the defaults. It is
   a submit button, so the page needs no JavaScript for it; the server re-renders it the other way. */
button.switch {{ display: inline-flex; align-items: center; gap: .5rem; }}
button.switch .track {{ position: relative; width: 2.1rem; height: 1.15rem; border-radius: 999px;
                        border: 1px solid var(--line); background: #fff; flex: none; }}
button.switch .knob {{ position: absolute; top: 1px; left: 1px; width: .95rem; height: .95rem;
                       border-radius: 50%; background: var(--muted); }}
button.switch.on .track {{ background: var(--made); border-color: var(--made); }}
button.switch.on .knob {{ left: auto; right: 1px; background: #fff; }}
/* the two input boxes side by side on a wide screen, one above the other on a narrow one */
.boxes {{ display: flex; gap: 1.2rem; flex-wrap: wrap; align-items: start; }}
.box {{ flex: 1 1 22rem; min-width: 0; display: flex; flex-direction: column; }}
/* the examples above one box may wrap where the other's do not, so the line keeps room for two either
   way and both text areas start at the same height */
.box .examples {{ min-height: 3.1em; }}
.backtop {{ justify-content: flex-end; margin: 0 0 -.4rem; }}
details {{ margin-top: 1.1rem; border-top: 1px solid var(--line); padding-top: .8rem; }}
summary {{ cursor: pointer; font-weight: 600; }}
.row {{ margin: .55rem 0; }}
.note {{ background: var(--panel); border: 1px solid var(--line); padding: .7rem .9rem; border-radius: 8px; }}
.result {{ margin-top: 1.4rem; border-top: 1px solid var(--line); padding-top: .4rem; }}
.scroll {{ overflow-x: auto; }} .legend svg {{ max-width: 100%; height: auto; }}
table {{ border-collapse: collapse; width: 100%; font-size: .9rem; margin-top: .5rem; }}
th, td {{ text-align: left; padding: .38rem .5rem; border-bottom: 1px solid var(--line); vertical-align: top; }}
th {{ color: var(--muted); font-weight: 600; font-size: .78rem; text-transform: uppercase; letter-spacing: .03em; }}
.nowrap {{ white-space: nowrap; }}
.made {{ color: var(--made); font-weight: 600; }} .took {{ color: var(--took); font-weight: 600; }}
.pill {{ display: inline-block; font-size: .75rem; padding: .05rem .45rem; margin: .1rem .15rem .1rem 0;
        border-radius: 999px; border: 1px solid var(--line); color: var(--muted); }}
.sources {{ font-size: .9rem; }} .sources li {{ margin-bottom: .3rem; }}
pre {{ background: var(--panel); border: 1px solid var(--line); border-radius: 8px; padding: .6rem .8rem;
      white-space: pre-wrap; }}
code {{ font-size: .9em; }} dd code {{ overflow-wrap: anywhere; }}
dt {{ font-weight: 600; margin-top: .8rem; }} dd {{ margin-left: 1rem; }}
.decisions li, .toc li {{ margin-bottom: .3rem; }}
a {{ color: var(--made); }}
form.inline {{ display: inline-flex; gap: .4rem; align-items: center; margin: 0; }}
.outputs {{ margin: .6rem 0 1rem; }}
.outputs details.report {{ margin: 0; border: 0; padding: 0; }}
.outputs details.report[open] {{ flex-basis: 100%; }}
summary.btn {{ list-style: none; display: inline-block; font-weight: normal; }}
summary.btn::-webkit-details-marker {{ display: none; }}
.report pre {{ max-height: 26rem; overflow: auto; font-size: .82rem; margin-top: .6rem; }}
progress {{ width: 100%; height: 10px; accent-color: var(--made); }}
/* the phase choice under the boxes: three radio buttons on one line */
.phase {{ display: flex; gap: 1.1rem; flex-wrap: wrap; align-items: center; margin-top: .8rem; }}
.phase label {{ display: inline-flex; gap: .35rem; align-items: center; }}
.na {{ color: var(--muted); font-style: italic; }}
td.num {{ text-align: right; font-variant-numeric: tabular-nums; }}
.figure svg {{ max-width: 100%; height: auto; display: block; }}
"""
