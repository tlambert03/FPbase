from __future__ import annotations

import io
import re
from xml.etree import ElementTree

import pytest
from PIL import Image

from proteins.util.spectra_plot import _shape_line, spectra_fig, spectra_figure, to_rgba


class FakeSpectrum:
    def __init__(self, lo: int = 300, hi: int = 800, color: str = "#ff4600") -> None:
        self.data = [(w, max(0.0, 1 - abs(w - 500) / 100)) for w in range(lo, hi + 1)]
        self.min_wave, self.max_wave = lo, hi
        self._color = color

    def color(self) -> str:
        return self._color


def _render(fmt: str, **kwargs) -> bytes:
    out = spectra_fig([FakeSpectrum()], fmt, **kwargs)
    assert out is not None
    return out.getvalue()


@pytest.mark.parametrize(
    ("fmt", "size"),
    [("png", (840, 210)), ("jpg", (840, 210)), ("tif", (840, 210)), ("tiff", (840, 210))],
)
def test_raster_formats(fmt: str, size: tuple[int, int]) -> None:
    image = Image.open(io.BytesIO(_render(fmt)))
    assert image.size == size
    assert image.mode == ("RGB" if fmt == "jpg" else "RGBA")


def test_twitter_png_is_opaque_and_twice_as_tall() -> None:
    image = Image.open(io.BytesIO(_render("png", twitter=1, title="EGFP", info="Ex/Em")))
    assert image.size == (840, 420)
    assert image.getchannel("A").getextrema() == (255, 255)


def test_svg_is_valid() -> None:
    root = ElementTree.fromstring(_render("svg", title="mCherry", grid=True, ylabels=True))
    assert root.get("width") == "864pt"
    assert root.get("height") == "216pt"


def test_pdf_is_well_formed() -> None:
    pdf = _render("pdf", title="mCherry")
    assert pdf.startswith(b"%PDF-1.4") and pdf.rstrip().endswith(b"%%EOF")
    xref = int(re.search(rb"startxref\n(\d+)", pdf).group(1))  # pyright: ignore
    assert pdf[xref:].startswith(b"xref")
    for n, offset in enumerate(re.findall(rb"(\d{10}) 00000 n", pdf), 1):
        assert pdf[int(offset) :].startswith(f"{n} 0 obj".encode())


def test_unsupported_format() -> None:
    with pytest.raises(ValueError, match="Unsupported image format"):
        _render("gif")


def test_no_spectra() -> None:
    assert spectra_fig([], "svg") is None


def test_text_layout_matches_matplotlib() -> None:
    # glyph positions (font units at size 100) from matplotlib 3.11's HarfBuzz layout
    *_, glyphs = _shape_line("office", 100)
    assert [round(x, 3) for _, x in glyphs] == [0, 61.188, 157.875, 212.859]  # ffi ligature
    *_, glyphs = _shape_line("AV", 100)
    assert round(glyphs[1][1], 4) == 62.0156  # kerned
    assert _shape_line("300", 100)[0] == 190.875
    assert _shape_line("mCherry", 100)[0] == 431.796875


def test_svg_text_positions_match_matplotlib() -> None:
    # anchors matplotlib 3.11 wrote for the same figure ("translate(x y) scale(...)")
    svg = _render("svg", title="mCherry", info="Ex/Em λ: 587/610\nEC: 72000  QY: 0.22")
    anchors = re.findall(r'transform="translate\(([\d.]+) ([\d.]+)\) scale', svg.decode())
    assert ("3.41625", "209.817656") in anchors  # "300" tick label
    assert ("16.31232", "25.31078") in anchors  # title
    assert ("16.31232", "46.698791") in anchors  # first info line


def test_ticks() -> None:
    fig = spectra_figure([FakeSpectrum()], ylabels=True, grid=True)
    labels = [t.text for t in fig.texts]
    assert labels[:11] == [str(x) for x in range(300, 801, 50)]
    assert labels[11:] == ["0.0", "0.2", "0.4", "0.6", "0.8", "1.0"]
    minor_ticks = [s for s in fig.shapes if s.width == 0.6]
    assert len(minor_ticks) == 40 + 5  # every 10 nm and 0.1, except at labeled ticks


@pytest.mark.parametrize(
    ("color", "rgba"),
    [
        ("#ff4600", (1.0, 70 / 255, 0.0, 1.0)),
        ("gray", (128 / 255, 128 / 255, 128 / 255, 1.0)),
        ("r", (1.0, 0.0, 0.0, 1.0)),
        ("0.5", (0.5, 0.5, 0.5, 1.0)),
    ],
)
def test_to_rgba(color: str, rgba: tuple[float, ...]) -> None:
    assert to_rgba(color) == pytest.approx(rgba)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"xlim": [0, 10**30]},  # would allocate ticks without bound
        {"xlim": [0, 10**5], "grid": True},  # >1000 minor ticks, matplotlib's limit too
        {"xlim": [0, float("inf")]},
        {"alpha": 5},
        {"alpha": "-1"},
        {"linewidth": -5},
    ],
)
def test_rejects_bad_options(kwargs: dict) -> None:
    with pytest.raises((ValueError, OverflowError)):
        _render("svg", **kwargs)


@pytest.mark.parametrize("fmt", ["svg", "png", "pdf"])
def test_accepts_edge_options(fmt: str) -> None:
    _render(fmt, xlim=[500, 500])  # widened, as matplotlib does
    _render(fmt, xlim=[0, 10**5], xlabels=False)  # no ticks to draw
    _render(fmt, linewidth=10**9, fill=False)  # clamped for rasters
