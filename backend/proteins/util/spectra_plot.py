"""Static spectra images (SVG, PNG, JPEG, TIFF, PDF) without matplotlib.

Reproduces the figures FPbase used to draw with matplotlib: same layout, ticks, colors
and DejaVu Sans text, so web workers don't need to load matplotlib.

- SVG and PDF draw text as glyph outlines placed exactly as matplotlib's (unhinted)
  layout places them, from a table extracted from matplotlib by
  `scripts/build_spectra_font.py`.
- Raster formats are drawn with Pillow, following matplotlib's Agg backend: shapes are
  antialiased (by supersampling), straight lines are snapped to pixels, and text is
  rendered hinted by FreeType at the hinted positions Agg uses.
"""

from __future__ import annotations

import io
import itertools
import json
import math
import re
import zlib
from dataclasses import dataclass, field
from functools import cache, partial
from pathlib import Path
from typing import TYPE_CHECKING, Any

from PIL import Image, ImageColor, ImageDraw, ImageFont

if TYPE_CHECKING:
    from collections.abc import Callable, Sequence

    from proteins.models import Spectrum

RGBA = tuple[float, float, float, float]
Point = tuple[float, float]
Box = tuple[float, float, float, float]  # x, y, width, height

FONT_DIR = Path(__file__).parent / "spectra_font"
DPI = 70  # raster resolution; figure sizes are in inches of 72 points
SUPERSAMPLE = 4  # raster shapes are drawn this much larger, then downsampled
TICK_GRAY = (0.2, 0.2, 0.2, 1.0)
SPINE_GRAY = (0.5, 0.5, 0.5, 1.0)
GRID = (0.5, 0.5, 0.5, 0.15)
BLACK = (0.0, 0.0, 0.0, 1.0)
FADED_BLACK = (0.0, 0.0, 0.0, 0.5)
WHITE = (1.0, 1.0, 1.0, 1.0)
TICK_LABEL_SIZE = 10
TICK_PAD = 3.5
MINOR_TICK_LENGTH, MINOR_TICK_WIDTH = 2.0, 0.6
FORMATS = {
    **{"svg": "svg", "pdf": "pdf", "png": "PNG", "jpg": "JPEG", "jpeg": "JPEG"},
    **{"tif": "TIFF", "tiff": "TIFF"},
}
MAX_TICKS = 1000  # per tick spacing, as matplotlib's `Locator.MAXTICKS`
# matplotlib's single-letter colors, which CSS doesn't have
BASE_COLORS = {"b": "#0000ff", "g": "#008000", "r": "#ff0000", "c": "#00bfbf"}
BASE_COLORS |= {"m": "#bf00bf", "y": "#bfbf00", "k": "#000000", "w": "#ffffff"}


def to_rgba(color: str, alpha: float = 1.0) -> RGBA:
    """Parse a hex, CSS or matplotlib-style color name, or a gray level like '0.5'."""
    color = BASE_COLORS.get(color, color)
    try:
        gray = float(color)
    except ValueError:
        r, g, b, *a = ImageColor.getrgb(color)
        return (r / 255, g / 255, b / 255, alpha * (a[0] / 255 if a else 1))
    return (gray, gray, gray, alpha)


# ------------------------------------------------------------------ figure model


@dataclass
class Shape:
    """A path in figure points (origin top left): filled and/or stroked."""

    points: list[Point]
    fill: RGBA | None = None
    stroke: RGBA | None = None
    width: float = 1.0
    cap: str = "butt"  # or "square" (matplotlib's "projecting")
    closed: bool = False
    clip: Box | None = None
    css_class: str = ""  # for the SVG


@dataclass
class Text:
    """Text anchored at `x, y` (points), aligned like matplotlib's `ha` and `va`."""

    text: str
    x: float
    y: float
    size: float
    color: RGBA
    ha: str = "left"  # left | center | right
    va: str = "top"  # top | center_baseline


@dataclass
class Figure:
    width: float  # points
    height: float
    background: RGBA | None
    shapes: list[Shape] = field(default_factory=list)
    texts: list[Text] = field(default_factory=list)


def _multiples(step: float, lo: float, hi: float) -> list[float]:
    """Multiples of `step` within [lo, hi], like matplotlib's `MultipleLocator`."""
    first, last = math.ceil(lo / step - 1e-9), math.floor(hi / step + 1e-9)
    if last - first >= MAX_TICKS:
        raise ValueError(f"x range too wide: would draw {last - first + 1} ticks")
    return [round(k * step, 10) for k in range(first, last + 1)]


def _auto_step(lo: float, hi: float, length: float) -> float:
    """The tick step `MaxNLocator` picks for a y axis `length` points tall."""
    nbins = min(max(int(length // (TICK_LABEL_SIZE * 2)), 1), 9)
    raw = (hi - lo) / nbins
    scale = 10 ** math.floor(math.log10(raw))
    return next(s * scale for s in (1, 2, 2.5, 5, 10) if s * scale >= raw)


def _tick_label(value: float, step: float) -> str:
    decimals = 0
    while abs(step * 10**decimals - round(step * 10**decimals)) > 1e-9:
        decimals += 1
    return f"{value:.{decimals}f}"


def spectra_figure(
    spectra: Sequence[Spectrum],
    xlabels: bool = True,
    ylabels: bool = False,
    xlim: Sequence[float] | None = None,
    fill: bool = True,
    transparent: bool = True,
    grid: bool = False,
    title: str | bool = False,
    info: str | None = None,
    figsize: tuple[float, float] = (12, 3),
    alpha: float | str | None = None,
    color: str | None = None,
    twitter: bool | int = False,
    linewidth: float | str | None = None,
    **_ignored: Any,
) -> Figure:
    """Lay out a spectra figure, as the former matplotlib `spectra_fig` did."""
    if twitter:
        transparent, figsize, xlim, xlabels = False, (12, 6), (400, 760), True
    if not xlim:
        xlim = (min(s.min_wave for s in spectra), max(s.max_wave for s in spectra))
    x0, x1 = float(xlim[0]), float(xlim[1])
    if x0 == x1:  # widen, as matplotlib's `nonsingular` does
        x0, x1 = (-0.05, 0.05) if x0 == 0 else (x0 - 0.05 * abs(x0), x1 + 0.05 * abs(x1))
    alpha = float(alpha) if alpha else None
    if alpha is not None and not 0 <= alpha <= 1:
        raise ValueError(f"alpha must be between 0 and 1, not {alpha}")
    linewidth = None if linewidth is None else float(linewidth)
    if linewidth is not None and linewidth < 0:
        raise ValueError(f"linewidth must not be negative, not {linewidth}")
    y0, y1 = (0, 1.07) if twitter else (-0.005, 1.025)

    W, H = figsize[0] * 72, figsize[1] * 72
    pos = [0, 0.017, 0.97, 0.98]  # axes left, bottom, width, height (fractions)
    if xlabels:
        pos[0], pos[1], pos[3] = 0.02, 0.08, pos[3] - 0.065
    if ylabels:
        pos[0], pos[2], pos[3] = 0.025, 0.96, pos[3] - 0.01
    else:
        pos[0] = 0.015
    left, width, height = pos[0] * W, pos[2] * W, pos[3] * H
    top = H - (pos[1] + pos[3]) * H
    bottom = top + height
    clip = (left, top, width, height)

    def X(x: float) -> float:
        return left + (x - x0) / (x1 - x0) * width

    def Y(y: float) -> float:
        return bottom - (y - y0) / (y1 - y0) * height

    fig = Figure(W, H, None if transparent else WHITE)
    shapes, texts = fig.shapes, fig.texts

    # matplotlib hides the major tick marks; minor ticks only exist with a grid
    ystep = _auto_step(y0, y1, height)
    xmajor = _multiples(50, x0, x1) if xlabels else []
    ymajor = _multiples(ystep, y0, y1) if ylabels else []
    if grid:  # drawn below everything else
        xminor = [x for x in _multiples(10, x0, x1) if x not in xmajor]
        yminor = [y for y in _multiples(0.1, y0, y1) if y not in ymajor]
        for x in xmajor + xminor:
            pts = [(X(x), bottom), (X(x), top)]
            shapes.append(Shape(pts, stroke=GRID, width=0.8, cap="square", clip=clip))
        for x in xminor:
            pts = [(X(x), bottom), (X(x), bottom + MINOR_TICK_LENGTH)]
            shapes.append(Shape(pts, stroke=BLACK, width=MINOR_TICK_WIDTH))
        for y in ymajor + yminor:
            pts = [(left, Y(y)), (left + width, Y(y))]
            shapes.append(Shape(pts, stroke=GRID, width=0.8, cap="square", clip=clip))
        for y in yminor:
            pts = [(left, Y(y)), (left - MINOR_TICK_LENGTH, Y(y))]
            shapes.append(Shape(pts, stroke=BLACK, width=MINOR_TICK_WIDTH))

    for spec in spectra:
        pts = [(X(x), Y(y)) for x, y in spec.data]
        if fill:
            face = to_rgba(color or spec.color(), alpha or 0.5)
            edge = 1.0 if linewidth is None else linewidth
            base = Y(0)
            outline = [(pts[0][0], base), *pts, *((x, base) for x, _ in reversed(pts))]
            edge_color = face if edge else None
            shapes.append(
                Shape(
                    outline, face, edge_color, edge, closed=True, clip=clip, css_class="spectrum"
                )
            )
        else:
            stroke = to_rgba(spec.color(), alpha or 1)
            lw = 1.5 if linewidth is None else linewidth
            shapes.append(
                Shape(pts, stroke=stroke, width=lw, cap="square", clip=clip, css_class="spectrum")
            )

    if xlabels:
        shapes.append(
            Shape([(left, bottom), (left + width, bottom)], None, SPINE_GRAY, 0.4, "square")
        )
        for x in xmajor:
            label = _tick_label(x, 50)
            texts.append(
                Text(label, X(x), bottom + TICK_PAD, TICK_LABEL_SIZE, TICK_GRAY, "center")
            )
    if ylabels:
        shapes.append(Shape([(left, bottom), (left, top)], None, SPINE_GRAY, 0.4, "square"))
        for y in ymajor:
            label, x = _tick_label(y, ystep), left - TICK_PAD
            texts.append(
                Text(label, x, Y(y), TICK_LABEL_SIZE, TICK_GRAY, "right", "center_baseline")
            )
    if title:
        texts.append(Text(str(title), X(x0 + 2), Y(0.97), 18, FADED_BLACK))
        if info:
            texts.append(Text(info, X(x0 + 2), Y(0.85), 14, FADED_BLACK))
    return fig


# ------------------------------------------------------------------ text layout


@cache
def _font_table() -> dict[str, Any]:
    return json.loads((FONT_DIR / "glyphs.json").read_text())


@cache
def _pil_font(size: float) -> ImageFont.FreeTypeFont:
    return ImageFont.truetype(FONT_DIR / "DejaVuSans.ttf", size)


# (width, ascent, descent, glyphs) of one line of text, in the caller's units
Metrics = tuple[float, float, float, Any]


def layout_text(
    text: Text,
    measure: Callable[[str], Metrics],
    min_ascent: float,
    min_descent: float,
    line_gap: float,
) -> list[tuple[float, float, Any]]:
    """Place each line of `text`, as `matplotlib.text.Text._get_layout` does.

    `measure` and the other metrics are in the units of `text.x, text.y`. Returns
    `(x, baseline, glyphs)` per line.
    """
    lines = text.text.split("\n")
    if len(lines) == 1:
        line_gap = 0
    placed, widths, cursor, first_ascent = [], [], 0.0, 0.0
    for i, line in enumerate(lines):
        width, ascent, descent, glyphs = measure(line) if line else (0, 0, 0, [])
        ascent = max(ascent, min_ascent) + line_gap / 2
        descent = max(descent, min_descent) + line_gap / 2
        if i == 0:
            first_ascent = ascent
        cursor += ascent
        placed.append((cursor, glyphs))
        widths.append(width)
        cursor += descent
    width = max(widths)
    dx = {"left": 0, "center": width / 2, "right": width}[text.ha]
    dy = {"top": 0, "center_baseline": first_ascent / 2}[text.va]
    return [(text.x - dx, text.y - dy + base, glyphs) for base, glyphs in placed]


def _tokenize(line: str) -> list[tuple[str, float]]:
    """Split `line` into glyphs (ligatures first, as HarfBuzz does), each with the pair
    kerning before it, in font units.
    """
    table = _font_table()
    glyphs, kerning = table["glyphs"], table["kerning"]
    tokens, i, prev = [], 0, None
    while i < len(line):
        n = next(n for n in (3, 2, 1) if line[i : i + n] in glyphs or n == 1)
        tok = line[i : i + n]
        tokens.append((tok, kerning.get(f"{prev}\0{tok}", 0) if prev else 0))
        prev, i = tok, i + n
    return tokens


def _shape_line(line: str, size: float) -> Metrics:
    """Measure and place glyph outlines (offsets in font units) as matplotlib's unhinted
    layout does. Metrics are in points.
    """
    table = _font_table()
    glyphs, x = [], 0.0
    ymin, ymax = math.inf, -math.inf
    for tok, kern in _tokenize(line):
        advance, lo, hi, path = table["glyphs"].get(tok, table["notdef"])
        glyphs.append((path, x + kern))
        ymin, ymax = min(ymin, lo), max(ymax, hi)
        x += kern + advance
    scale = size / table["font_scale"]
    return x * scale, ymax * scale, -ymin * scale, glyphs


@dataclass
class GlyphRun:
    """One line of glyph outlines, in font units, at a baseline in points."""

    glyphs: list[tuple[str, float]]  # (path, x offset)
    x: float
    y: float
    size: float
    color: RGBA

    @property
    def scale(self) -> float:
        return self.size / _font_table()["font_scale"]


def _layout(text: Text, measure: Callable[[str], Metrics]) -> list[tuple[float, float, Any]]:
    """`layout_text` with DejaVu Sans' minimum line metrics, as matplotlib uses."""
    table, size = _font_table(), text.size
    ascent, descent, gap = (table[k] * size for k in ("ascent", "descent", "line_gap"))
    return layout_text(text, measure, ascent, descent, gap)


def _glyph_runs(fig: Figure) -> list[GlyphRun]:
    return [
        GlyphRun(glyphs, x, y, t.size, t.color)
        for t in fig.texts
        for x, y, glyphs in _layout(t, partial(_shape_line, size=t.size))
    ]


# ------------------------------------------------------------------ SVG


def _num(v: float) -> str:
    return f"{v:.6f}".rstrip("0").rstrip(".")


def _hex(c: RGBA) -> str:
    return "#{:02x}{:02x}{:02x}".format(*(round(v * 255) for v in c[:3]))


def render_svg(fig: Figure) -> bytes:
    w, h = _num(fig.width), _num(fig.height)
    out = [
        '<?xml version="1.0" encoding="utf-8" standalone="no"?>\n',
        '<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink"',
        f' width="{w}pt" height="{h}pt" viewBox="0 0 {w} {h}" version="1.1">\n',
        "<defs><style>*{stroke-linejoin:round;stroke-linecap:butt}</style>\n",
    ]
    clips: dict[Box, str] = {}
    for shape in fig.shapes:
        if shape.clip and shape.clip not in clips:
            clips[shape.clip] = f"clip{len(clips)}"
            x, y, cw, ch = map(_num, shape.clip)
            out.append(
                f'<clipPath id="{clips[shape.clip]}">'
                f'<rect x="{x}" y="{y}" width="{cw}" height="{ch}"/></clipPath>\n'
            )
    runs = _glyph_runs(fig)
    glyph_ids: dict[str, str] = {}
    for run in runs:
        for path, _ in run.glyphs:
            if path not in glyph_ids:
                glyph_ids[path] = gid = f"g{len(glyph_ids)}"
                out.append(f'<path id="{gid}" d="{path}" transform="scale(0.015625)"/>\n')
    out.append("</defs>\n")
    if fig.background:
        out.append(f'<rect width="{w}" height="{h}" style="fill:{_hex(fig.background)}"/>\n')
    for s in fig.shapes:
        d = (
            "M"
            + " L".join(f"{_num(x)} {_num(y)}" for x, y in s.points)
            + (" Z" if s.closed else "")
        )
        style = [f"fill:{_hex(s.fill)}" if s.fill else "fill:none"]
        if s.fill and s.fill[3] < 1:
            style.append(f"fill-opacity:{_num(s.fill[3])}")
        if s.stroke and s.width > 0:
            style.append(f"stroke:{_hex(s.stroke)};stroke-width:{_num(s.width)}")
            if s.stroke[3] < 1:
                style.append(f"stroke-opacity:{_num(s.stroke[3])}")
            if s.cap != "butt":
                style.append(f"stroke-linecap:{s.cap}")
        attrs = f' class="{s.css_class}"' if s.css_class else ""
        attrs += f' clip-path="url(#{clips[s.clip]})"' if s.clip else ""
        out.append(f'<path d="{d}"{attrs} style="{";".join(style)}"/>\n')
    for run in runs:
        style = f"fill:{_hex(run.color)}" + (
            f";opacity:{_num(run.color[3])}" if run.color[3] < 1 else ""
        )
        k = _num(run.scale)
        transform = f"translate({_num(run.x)} {_num(run.y)}) scale({k} -{k})"
        out.append(f'<g style="{style}" transform="{transform}">')
        for path, x in run.glyphs:
            offset = f' x="{_num(x)}"' if x else ""
            out.append(f'<use xlink:href="#{glyph_ids[path]}"{offset}/>')
        out.append("</g>\n")
    out.append("</svg>\n")
    return "".join(out).encode()


# ------------------------------------------------------------------ PDF

_PATH_TOKENS = re.compile(r"[MLQCZ]|-?\d+")


def _pdf_glyph_ops(path: str, run: GlyphRun, x0: float, page_height: float) -> str:
    """PDF path operators for a glyph (quadratic curves become cubic ones)."""
    k = run.scale / 64

    def pt(gx: float, gy: float) -> str:
        x, y = run.x + x0 * run.scale + gx * k, page_height - run.y + gy * k
        return f"{_num(x)} {_num(y)}"

    ops, cur, tokens, i = [], (0, 0), _PATH_TOKENS.findall(path), 0
    while i < len(tokens):
        cmd, nums = tokens[i], []
        i += 1
        while i < len(tokens) and tokens[i] not in "MLQCZ":
            nums.append(int(tokens[i]))
            i += 1
        pts = list(zip(nums[::2], nums[1::2], strict=True))
        if cmd == "M":
            ops.append(f"{pt(*pts[0])} m")
        elif cmd == "L":
            ops.append(f"{pt(*pts[0])} l")
        elif cmd == "Q":
            (cx, cy), end = pts
            c1 = (cur[0] + 2 / 3 * (cx - cur[0]), cur[1] + 2 / 3 * (cy - cur[1]))
            c2 = (end[0] + 2 / 3 * (cx - end[0]), end[1] + 2 / 3 * (cy - end[1]))
            ops.append(f"{pt(*c1)} {pt(*c2)} {pt(*end)} c")
        elif cmd == "C":
            ops.append(" ".join(pt(*p) for p in pts) + " c")
        elif cmd == "Z":
            ops.append("h")
        if pts:
            cur = pts[-1]
    return " ".join(ops)


def render_pdf(fig: Figure) -> bytes:
    """A one-page vector PDF (text as glyph outlines, so no fonts are embedded)."""
    H = fig.height
    ops: list[str] = []
    alphas: dict[float, str] = {}

    def alpha(a: float) -> str:
        alphas.setdefault(a, f"GS{len(alphas)}")
        return f"/{alphas[a]} gs"

    def rgb(c: RGBA) -> str:
        return " ".join(_num(v) for v in c[:3])

    def path(points: list[Point], closed: bool) -> str:
        cmds = [f"{_num(x)} {_num(H - y)} {'l' if i else 'm'}" for i, (x, y) in enumerate(points)]
        return " ".join(cmds) + (" h" if closed else "")

    if fig.background:
        ops.append(f"{rgb(fig.background)} rg 0 0 {_num(fig.width)} {_num(H)} re f")
    for s in fig.shapes:
        ops.append("q")
        if s.clip:
            x, y, w, h = s.clip
            ops.append(f"{_num(x)} {_num(H - y - h)} {_num(w)} {_num(h)} re W n")
        if s.fill:
            ops.append(f"{alpha(s.fill[3])} {rgb(s.fill)} rg {path(s.points, s.closed)} f")
        if s.stroke and s.width > 0:
            cap = {"butt": 0, "square": 2}[s.cap]
            ops.append(
                f"{alpha(s.stroke[3])} {rgb(s.stroke)} RG {_num(s.width)} w {cap} J 1 j "
                f"{path(s.points, s.closed)} S"
            )
        ops.append("Q")
    for run in _glyph_runs(fig):
        ops.append(f"q {alpha(run.color[3])} {rgb(run.color)} rg")
        ops += [_pdf_glyph_ops(p, run, x0, H) for p, x0 in run.glyphs]
        ops.append("f Q")

    content = zlib.compress("\n".join(ops).encode())
    gstates = " ".join(f"/{n} << /ca {_num(a)} /CA {_num(a)} >>" for a, n in alphas.items())
    objects = [
        b"<< /Type /Catalog /Pages 2 0 R >>",
        b"<< /Type /Pages /Kids [3 0 R] /Count 1 >>",
        (
            f"<< /Type /Page /Parent 2 0 R /MediaBox [0 0 {_num(fig.width)} {_num(H)}]"
            f" /Resources << /ExtGState << {gstates} >> >> /Contents 4 0 R >>"
        ).encode(),
        f"<< /Length {len(content)} /Filter /FlateDecode >>\nstream\n".encode()
        + content
        + b"\nendstream",
    ]
    pdf, offsets = bytearray(b"%PDF-1.4\n"), []
    for n, obj in enumerate(objects, 1):
        offsets.append(len(pdf))
        pdf += f"{n} 0 obj\n".encode() + obj + b"\nendobj\n"
    xref = len(pdf)
    pdf += f"xref\n0 {len(objects) + 1}\n0000000000 65535 f \n".encode()
    pdf += "".join(f"{o:010d} 00000 n \n" for o in offsets).encode()
    pdf += f"trailer\n<< /Size {len(objects) + 1} /Root 1 0 R >>\n".encode()
    pdf += f"startxref\n{xref}\n%%EOF\n".encode()
    return bytes(pdf)


# ------------------------------------------------------------------ raster


def _snap(points: list[Point], width: float) -> list[Point]:
    """Agg's `PathSnapper`: put horizontal/vertical paths on pixel centers or edges."""
    for (ax, ay), (bx, by) in itertools.pairwise(points):
        if abs(ax - bx) >= 1e-4 and abs(ay - by) >= 1e-4:
            return points
    offset = 0.5 if math.floor(width + 0.5) % 2 else 0.0
    return [(math.floor(x + 0.5) + offset, math.floor(y + 0.5) + offset) for x, y in points]


def _project_ends(points: list[Point], d: float) -> list[Point]:
    """Lengthen an open path by `d` at both ends (a square line cap)."""

    def push(p: Point, q: Point) -> Point:
        dx, dy = p[0] - q[0], p[1] - q[1]
        n = math.hypot(dx, dy) or 1
        return (p[0] + dx / n * d, p[1] + dy / n * d)

    return [push(points[0], points[1]), *points[1:-1], push(points[-1], points[-2])]


def _paint(canvas: Image.Image, mask: Image.Image, at: tuple[int, int], color: RGBA) -> None:
    rgb = tuple(round(v * 255) for v in color[:3])
    layer = Image.new("RGBA", mask.size, rgb)  # pyright: ignore[reportArgumentType]
    layer.putalpha(mask.point(lambda v: round(v * color[3])))
    canvas.alpha_composite(layer, at)


def _draw_shapes(canvas: Image.Image, fig: Figure, px: float) -> None:
    """Paint antialiased shapes, like Agg: each shape's coverage (here from a
    `SUPERSAMPLE`x larger mask, box-downsampled) is composited onto the canvas.
    """
    ss = SUPERSAMPLE
    for s in fig.shapes:
        # pixel coordinates (y down), then supersampled with Pillow's pixel centers
        points = [(x * px, y * px) for x, y in s.points]
        # strokes wider than the canvas all look alike (and Pillow's cost grows with width)
        width = min(s.width * px, 2 * sum(canvas.size))
        if s.stroke and not s.fill:
            points = _snap(points, width)
        bounds = (0, 0, *canvas.size)
        if s.clip:  # Agg rounds the clip box to whole pixels
            x, y, w, h = (v * px for v in s.clip)
            bounds = (round(x), round(y), round(x + w), round(y + h))
        pad = math.ceil(width) + 2
        xs, ys = [p[0] for p in points], [p[1] for p in points]
        x0 = max(bounds[0], math.floor(min(xs)) - pad)
        y0 = max(bounds[1], math.floor(min(ys)) - pad)
        x1 = min(bounds[2], math.ceil(max(xs)) + pad)
        y1 = min(bounds[3], math.ceil(max(ys)) + pad)
        if x1 <= x0 or y1 <= y0:
            continue
        local = [((x - x0) * ss - 0.5, (y - y0) * ss - 0.5) for x, y in points]
        mask_size = ((x1 - x0) * ss, (y1 - y0) * ss)
        if s.fill:
            # Agg draws matplotlib's fill_between faces (not their edges) half a pixel
            # right and down of the path; measured, as the source doesn't make it obvious
            face = [(x + ss / 2, y + ss / 2) for x, y in local]
            mask = Image.new("L", mask_size, 0)
            ImageDraw.Draw(mask).polygon(face, fill=255)
            _paint(canvas, mask.reduce(ss), (x0, y0), s.fill)
        if s.stroke and s.width > 0:
            mask = Image.new("L", mask_size, 0)
            line_width = max(1, round(width * ss))
            if s.closed:
                local = [*local, local[0]]
            elif s.cap == "square":
                local = _project_ends(local, line_width / 2)
            ImageDraw.Draw(mask).line(local, fill=255, width=line_width, joint="curve")
            _paint(canvas, mask.reduce(ss), (x0, y0), s.stroke)


def _measure_hinted(line: str, font: ImageFont.FreeTypeFont) -> Metrics:
    """Measure and place glyphs (offsets in pixels) as Agg does: by their hinted, whole
    pixel advances, plus unhinted kerning.
    """
    scale = font.size / _font_table()["font_scale"]
    glyphs, x = [], 0.0
    for tok, kern in _tokenize(line):
        glyphs.append((tok, x + kern * scale))
        x += kern * scale + font.getlength(tok)
    _, top, _, bottom = font.getbbox(line, anchor="ls")
    return x, -top, bottom, glyphs


def render_raster(fig: Figure, dpi: float = DPI) -> Image.Image:
    """Render to an RGBA image the way matplotlib's Agg backend would."""
    px = dpi / 72
    size = (round(fig.width * px), round(fig.height * px))
    canvas = Image.new("RGBA", size, (255, 255, 255, 0))
    if fig.background:
        background = tuple(round(v * 255) for v in fig.background)
        canvas.paste(background, (0, 0, *size))  # pyright: ignore[reportArgumentType]
    _draw_shapes(canvas, fig, px)

    masks: dict[RGBA, Image.Image] = {}  # one per color: texts don't overlap
    for t in fig.texts:
        font = _pil_font(t.size * px)
        in_pixels = Text(t.text, t.x * px, t.y * px, t.size * px, t.color, t.ha, t.va)
        if t.color not in masks:
            masks[t.color] = Image.new("L", size, 0)
        draw = ImageDraw.Draw(masks[t.color])
        for x, y, glyphs in _layout(in_pixels, partial(_measure_hinted, font=font)):
            for tok, dx in glyphs:
                draw.text((x + dx, y), tok, fill=255, font=font, anchor="ls")
    for color, mask in masks.items():
        _paint(canvas, mask, (0, 0), color)
    return canvas


# ------------------------------------------------------------------ entry point


def spectra_fig(
    spectra: Sequence[Spectrum],
    format: str = "svg",
    output: io.BytesIO | None = None,
    **kwargs: Any,
) -> io.BytesIO | None:
    """Draw `spectra` into `output` (a new `BytesIO` by default) as `format`."""
    if not spectra:
        return None
    if (fmt := FORMATS.get(format.lower())) is None:
        raise ValueError(f"Unsupported image format: {format!r}")
    fig = spectra_figure(spectra, **kwargs)
    output = output or io.BytesIO()
    if fmt == "svg":
        output.write(render_svg(fig))
    elif fmt == "pdf":
        output.write(render_pdf(fig))
    else:
        image = render_raster(fig)
        if fmt == "JPEG":  # no alpha channel: flatten onto white, as matplotlib does
            flat = Image.new("RGB", image.size, "white")
            flat.paste(image, mask=image.getchannel("A"))
            image = flat
        image.save(output, format=fmt, dpi=(DPI, DPI))
    return output
