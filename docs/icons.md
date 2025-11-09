# FPbase Icon System

FPbase uses an abstracted icon system that renders inline SVG icons from semantic icon names. This eliminates the need for external icon font CDNs while providing flexibility to switch between icon libraries.

## Overview

The icon system consists of three main components:

1. **Extraction Scripts** - Convert icon libraries (FontAwesome, Lucide, etc.) to a standard JSON format
2. **Icon Data** - JSON files containing SVG path data for each icon
3. **Template Tag** - Django template tag that renders inline SVG from icon data

## Using Icons in Templates

Use the `{% icon %}` template tag with semantic icon names:

```django
{% load fpbase_tags %}

{# Basic icon #}
{% icon "external-link" %}

{# Icon with CSS class #}
{% icon "info" class_="mr-2" %}

{# Icon with inline style #}
{% icon "warning" style="color: red;" %}

{# Icon with multiple attributes #}
{% icon "download" class_="icon-lg" aria_hidden="true" title="Download" %}
```

## Available Icons

The following semantic icon names are available:

- **UI & Navigation**: info, warning, alert, help, question, close, remove, menu, grid, search, filter, view, settings, edit, delete, trash, undo, check, success, selected, unselected, heart
- **Actions**: add, add-item, download, upload, share, share-square, link, external-link, exchange
- **Content**: book, collection, quote, photo, chart, table, flag, flag-outline
- **Time & Status**: clock, spinner, lightbulb, sun
- **Communication**: email
- **Tools**: wrench, keyboard
- **Social/External**: google, twitter, orcid

## Switching Icon Libraries

FPbase currently supports two icon libraries:

1. **FontAwesome** (default) - Filled icons, supports brand icons
2. **Lucide** - Outline/stroke icons, no brand icons

To switch libraries, set the `FPBASE_ICON_LIBRARY` environment variable:

```bash
export FPBASE_ICON_LIBRARY=lucide
```

Or update `backend/config/settings/base.py`:

```python
FPBASE_ICON_LIBRARY = "lucide"
```

## Adding a New Icon Library

To add support for a new icon library:

### 1. Install the library

```bash
pnpm add -D -w <library-name>
```

### 2. Create an extractor

Create a new extractor in `scripts/icon_extractors/<library>.py`:

```python
from __future__ import annotations

from pathlib import Path
from . import IconData, IconExtractor

class MyLibraryExtractor(IconExtractor):
    def __init__(self, base_path: Path | None = None):
        if base_path is None:
            base_path = Path(__file__).parent.parent.parent / "node_modules" / "<library>"
        self.base_path = base_path
        # Initialize paths...

    def extract_icon(self, icon_name: str) -> IconData | None:
        # Extract icon SVG data and return as IconData
        pass

    def get_available_icons(self) -> list[str]:
        # Return list of all available icon names
        pass
```

### 3. Update the extraction script

Add your library mapping to `scripts/extract_icons.py`:

```python
MY_LIBRARY_ICON_MAP = {
    "info": "information-icon",
    "warning": "alert-triangle",
    # ... map all FPbase icon names to your library's names
}

def extract_my_library_icons(output_dir: Path) -> None:
    extractor = MyLibraryExtractor()
    # ... extraction logic
```

Update the CLI to support your library:

```python
parser.add_argument(
    "library",
    choices=["fontawesome", "lucide", "mylibrary", "all"],
    # ...
)
```

### 4. Extract the icons

```bash
python scripts/extract_icons.py mylibrary
```

This will create `backend/fpbase/static/icons/mylibrary.json`.

### 5. Use the new library

Set the environment variable or update settings:

```bash
export FPBASE_ICON_LIBRARY=mylibrary
```

## Icon Data Format

All icon libraries are converted to a standard JSON format:

```json
{
  "external-link": {
    "viewBox": "0 0 24 24",
    "paths": [
      {
        "d": "M15 3h6v6...",
        "fill": "currentColor"
      }
    ]
  }
}
```

### IconPath Attributes

The `paths` array contains objects with the following optional attributes:

- `d` (required) - SVG path data
- `fill` - Fill color (usually "currentColor" or "none")
- `fillRule` - Fill rule (e.g., "evenodd")
- `clipRule` - Clip rule
- `stroke` - Stroke color (for stroke-based icons like Lucide)
- `strokeWidth` - Stroke width
- `strokeLinecap` - Stroke linecap style
- `strokeLinejoin` - Stroke linejoin style

## Adding New Icons

To add a new semantic icon to FPbase:

### 1. Add to AVAILABLE_ICONS

Update `backend/fpbase/templatetags/fpbase_tags.py`:

```python
AVAILABLE_ICONS = {
    # ... existing icons
    "my-new-icon",
}
```

### 2. Map to library icons

Update the icon maps in `scripts/extract_icons.py`:

```python
FONTAWESOME_ICON_MAP = {
    # ... existing mappings
    "my-new-icon": ("fas", "my-fa-icon"),
}

LUCIDE_ICON_MAP = {
    # ... existing mappings
    "my-new-icon": "my-lucide-icon",
}
```

### 3. Re-extract icons

```bash
python scripts/extract_icons.py all
```

### 4. Update tests

The tests will automatically include your new icon (via `AVAILABLE_ICONS`), but you may want to add specific tests.

## Bundle Size Comparison

- **Before**: ~900KB FontAwesome CDN (30,000+ icons)
- **After**: ~15KB inline SVG (49 icons)
- **Savings**: ~98% reduction in icon-related page weight

## Architecture Benefits

1. **Zero External Dependencies** - No CDN required, works offline
2. **Tree-Shakeable** - Only includes icons actually used
3. **Library Agnostic** - Easy to switch between icon libraries
4. **Type Safe** - Icon names validated at template render time
5. **XSS Protection** - All user input escaped
6. **Future Proof** - Adding new libraries is straightforward
