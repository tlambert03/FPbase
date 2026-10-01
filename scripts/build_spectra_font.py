"""Build `backend/proteins/util/spectra_font/` for the spectra image renderer.

SVG and PDF output draw text as glyph outlines, placed as matplotlib's (HarfBuzz)
layout places them, so `glyphs.json` holds DejaVu Sans outlines, advances, pair
kerning and ligatures, all extracted from matplotlib (which is otherwise not a
dependency). Raster output renders text with FreeType via Pillow, from a subset of
matplotlib's copy of `DejaVuSans.ttf`.

    uv run --with matplotlib --with fonttools python scripts/build_spectra_font.py
"""

from __future__ import annotations

import json
import string
from itertools import product
from pathlib import Path

from fontTools import subset
from matplotlib import _text_helpers
from matplotlib.font_manager import FontProperties, findfont, get_font
from matplotlib.ft2font import LoadFlags
from matplotlib.path import Path as MplPath

OUT = Path(__file__).parents[1] / "backend/proteins/util/spectra_font"
# matplotlib path code: (SVG command, number of vertices)
PATH_COMMANDS = {
    MplPath.MOVETO: ("M", 1),
    MplPath.LINETO: ("L", 1),
    MplPath.CURVE3: ("Q", 2),
    MplPath.CURVE4: ("C", 3),
    MplPath.CLOSEPOLY: ("Z", 1),
}
FONT_SCALE = 100  # matplotlib's TextToPath.FONT_SCALE; glyph paths are in 1/64 units

CHARS = sorted(
    set(map(chr, range(0x20, 0x7F)))  # ASCII
    | set(map(chr, range(0xA0, 0x180)))  # Latin-1 supplement, Latin extended-A
    | set(map(chr, range(0x391, 0x3CA)))  # Greek
    | set(map(chr, range(0x1EA0, 0x1EFA)))  # Vietnamese
    # punctuation, symbols, superscripts and subscripts
    | set(
        "\u2010\u2013\u2014\u2018\u2019\u201c\u201d\u2022\u2026\u2032\u2033\u00b1\u00d7\u2212\u2192\u2190\u2194\u2206\u00b5"
    )
    | set(
        "\u2070\u00b9\u00b2\u00b3\u2074\u2075\u2076\u2077\u2078\u2079\u2080\u2081\u2082\u2083\u2084\u2085\u2086\u2087\u2088\u2089"
    )
)

font = get_font(findfont(FontProperties(family=["sans-serif"])))
assert Path(font.fname).name == "DejaVuSans.ttf", font.fname
font.set_size(FONT_SCALE, 72)


def layout(s: str) -> list[tuple[int, float]]:
    return [(item.glyph_index, item.x) for item in _text_helpers.layout(s, font)]


def pen_end(s: str) -> float:
    font.set_text(s, 0, flags=LoadFlags.NO_HINTING)
    return font.get_width_height()[0] / 64


def glyph(index: int, text: str) -> list:
    font.set_text(text, 0, flags=LoadFlags.NO_HINTING)
    _, h = font.get_width_height()
    descent = font.get_descent()
    font.load_glyph(index, flags=LoadFlags.NO_HINTING)
    verts, codes = font.get_path()
    verts = (verts * 64).round().astype(int)  # back to FreeType's 26.6 units
    cmds, i = [], 0
    while i < len(codes):  # codes has one entry per vertex
        code = codes[i]
        cmd, n = PATH_COMMANDS[code]
        nums = [] if cmd == "Z" else verts[i : i + n].ravel().tolist()
        cmds.append(cmd + " ".join(map(str, nums)))
        i += n
    # [advance, ink ymin, ink ymax, path]; ink extents in FONT_SCALE units
    return [pen_end(text), -descent / 64, (h - descent) / 64, "".join(cmds)]


glyphs: dict[str, list] = {}
for c in CHARS:
    (index, _), *rest = layout(c)
    assert not rest
    if index:  # 0 is .notdef: the font lacks this character
        glyphs[c] = glyph(index, c)

# ligatures: runs of 2-3 letters that shape to a single glyph
ligatures = {}
letters = string.ascii_letters
for n in (2, 3):
    for chars in product(letters, repeat=n):
        s = "".join(chars)
        if len(items := layout(s)) == 1:
            ligatures[s] = glyph(items[0][0], s)
assert all(lig[:-1] in ligatures or len(lig) == 2 for lig in ligatures), ligatures
glyphs.update(ligatures)

kerning = {}
tokens = list(glyphs)
for a in tokens:
    for b in tokens:
        if a in ligatures and b in ligatures:
            continue
        items = layout(a + b)
        if len(items) == 2 and abs(kern := items[1][1] - glyphs[a][0]) > 1e-3:
            kerning[a + "\0" + b] = kern

os2 = font.get_sfnt_table("OS/2")
upem = font.get_sfnt_table("head")["unitsPerEm"]
data = {
    "font": "DejaVu Sans (glyph outlines from matplotlib's copy; Bitstream Vera license)",
    "font_scale": FONT_SCALE,
    # matplotlib's minimum line metrics, as fractions of the font size
    "ascent": os2["sTypoAscender"] / upem,
    "descent": -os2["sTypoDescender"] / upem,
    "line_gap": os2["sTypoLineGap"] / upem,
    "notdef": glyph(0, "￿"),
    "glyphs": glyphs,
    "kerning": kerning,
}
OUT.mkdir(exist_ok=True)
(OUT / "glyphs.json").write_text(
    json.dumps(data, ensure_ascii=False, separators=(",", ":")) + "\n"
)

# keep TrueType hinting and the legacy 'kern' table, which FreeType uses
options = subset.Options()
options.hinting = True
options.legacy_kern = True
options.layout_features = ["*"]
options.notdef_outline = True
ttf = subset.load_font(font.fname, options)
subsetter = subset.Subsetter(options)
subsetter.populate(unicodes=[ord(c) for c in CHARS])
subsetter.subset(ttf)
ttf.recalcTimestamp = False  # reproducible output
ttf.save(OUT / "DejaVuSans.ttf")
license_text = Path(font.fname).with_name("LICENSE_DEJAVU").read_text()
(OUT / "LICENSE_DEJAVU").write_text(
    "\n".join(line.rstrip() for line in license_text.splitlines()) + "\n"
)

print(f"{len(glyphs)} glyphs, {len(ligatures)} ligatures, {len(kerning)} kerning pairs")
for f in sorted(OUT.iterdir()):
    print(f"wrote {f} ({f.stat().st_size / 1024:.0f} KiB)")
