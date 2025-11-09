# extract_fontawesome_svgs.py
import argparse
import re
import shutil
import tempfile
import zipfile
from pathlib import Path

import requests

VERSION = "7.1.0"
URL = f"https://github.com/FortAwesome/Font-Awesome/releases/download/{VERSION}/fontawesome-free-{VERSION}-web.zip"

ICON_MAP = {
    # UI & Navigation
    "info": ("solid", "info-circle"),
    "info-i": ("solid", "info"),
    "warning": ("solid", "exclamation-circle"),
    "alert": ("solid", "exclamation-triangle"),
    "help": ("solid", "info-circle"),
    "question": ("solid", "question-circle"),
    "close": ("solid", "times"),
    "remove": ("solid", "times-circle"),
    "menu": ("solid", "list"),
    "grid": ("solid", "th"),
    "search": ("solid", "search"),
    "filter": ("solid", "filter"),
    "view": ("solid", "eye"),
    "settings": ("solid", "cog"),
    "edit": ("solid", "edit"),
    "delete": ("solid", "trash-alt"),
    "trash": ("solid", "trash"),
    "undo": ("solid", "undo"),
    "check": ("solid", "check"),
    "success": ("solid", "check-circle"),
    "selected": ("regular", "check-square"),
    "unselected": ("regular", "square"),
    "heart": ("solid", "heart"),
    "heart-outline": ("regular", "heart"),
    "add": ("solid", "plus"),
    "add-item": ("solid", "plus-circle"),
    "download": ("solid", "download"),
    "upload": ("solid", "upload"),
    "share": ("solid", "share"),
    "share-square": ("solid", "share-square"),
    "link": ("solid", "link"),
    "external-link": ("solid", "external-link-alt"),
    "exchange": ("solid", "exchange-alt"),
    "book": ("solid", "book"),
    "collection": ("solid", "book"),
    "quote": ("solid", "quote-left"),
    "photo": ("solid", "camera"),
    "chart": ("solid", "chart-area"),
    "table": ("solid", "table"),
    "flag": ("solid", "flag"),
    "flag-outline": ("regular", "flag"),
    "clock": ("solid", "clock"),
    "spinner": ("solid", "spinner"),
    "lightbulb": ("solid", "lightbulb"),
    "sun": ("solid", "sun"),
    "email": ("solid", "envelope"),
    "wrench": ("solid", "wrench"),
    "keyboard": ("regular", "keyboard"),
    "google": ("brands", "google"),
    "twitter": ("brands", "twitter"),
    "orcid": ("brands", "orcid"),
}


XML_COMMENT_RE = re.compile(r"<!--.*?-->", re.DOTALL)


def _copy_icon(src: Path, dest: Path) -> None:
    content = src.read_text(encoding="utf-8")
    content = XML_COMMENT_RE.sub("", content)
    dest.write_text(content, encoding="utf-8")


def extract_svgs(
    outdir: str | Path = "svgs", extract_dir_name: str = "fa-extract", keep_existing: bool = False
) -> None:
    """Download Font Awesome zip, extract into a named temp subfolder, and copy selected SVGs.

    Args:
        outdir: destination folder (created if missing) where the selected SVGs will be copied.
        extract_dir_name: the subfolder name under the temporary directory where the zip will be extracted.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        extract_path = tmpdir_path / extract_dir_name
        extract_path.mkdir(parents=True, exist_ok=True)
        zip_path = tmpdir_path / "fa.zip"

        print(f"Downloading {URL}...")
        with requests.get(URL, stream=True) as r:
            r.raise_for_status()
            with open(zip_path, "wb") as f:
                for chunk in r.iter_content(8192):
                    f.write(chunk)

        print("Extracting...")
        with zipfile.ZipFile(zip_path) as z:
            # Extract into a named subfolder so caller can control where files land
            z.extractall(extract_path)

        # Find top-level extracted Font-Awesome-* folder
        try:
            svg_path = next(extract_path.glob("fontawesome-*/svgs/"))
        except StopIteration:
            raise RuntimeError(
                "Could not find extracted Font Awesome folder. Check the structure of the zip file."
            ) from None

        outdir = Path(outdir)
        if outdir.exists() and not keep_existing:
            print("Removing existing SVGs...")
            shutil.rmtree(outdir, ignore_errors=True)
        outdir.mkdir(parents=True, exist_ok=True)

        for key, (style, glyph) in ICON_MAP.items():
            src = svg_path / style / f"{glyph}.svg"
            if src.exists():
                _copy_icon(src, outdir / f"{key}.svg")
                print(f"Copied {style}/{src.name} -> {key}.svg")
            else:
                raise ValueError(f"Missing: {glyph} ({style})")

    print("\n✅ Done")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download Font Awesome and extract selected SVGs.")
    parser.add_argument(
        "--outdir", "-o", default="backend/fpbase/static/icons", help="Destination folder for copied SVGs"
    )
    parser.add_argument(
        "--extract-dir-name",
        "-e",
        default="fa-extract",
        help="Name of temporary subfolder to extract the zip into",
    )
    parser.add_argument(
        "--keep-existing",
        action="store_true",
        help="If set, existing SVGs in the output directory will not be overwritten",
    )
    args = parser.parse_args()

    extract_svgs(outdir=args.outdir, extract_dir_name=args.extract_dir_name, keep_existing=args.keep_existing)
