# Legacy Static JavaScript

This directory contains legacy JavaScript that is **isolated from the main vite bundle** and loaded separately via Django static files.

## Contents

### `microscope.js` - Legacy Microscope Visualization (Isolated)

## Why is this isolated?

`microscope.js` is a complex 2100-line visualization using legacy libraries:

- **D3 v3.5.17** (modern codebase uses D3 v7.9.0)
- **NVD3 v1.8.6** (removed from modern codebase)
- **noUiSlider v10.1.0** (modern codebase bundles newer version)

Migrating this file to modern libraries would require 4-5 days of work due to:

- Custom brush/focus interactions tied to nvd3 API
- Dynamic scale switching (log/linear)
- Efficiency calculations that regenerate chart data
- Y-axis domain slider integration
- Mobile-specific handling
- 36+ chart method calls throughout the codebase

## How does it work?

### Build Process

1. Webpack `CopyPlugin` copies this file to `dist/js/microscope.js`
2. Django serves it as a static asset via `{% static 'js/microscope.js' %}`

### Page Loading

**Microscope detail pages** (`microscope_detail.html`, `microscope_embed.html`):

1. Load D3 v3, NVD3, noUiSlider from CDN (in `<head>`)
2. Load main webpack bundle (has D3 v7, modern code)
3. Load `microscope.js` script (uses CDN globals)

### Isolation Strategy

- ✅ Main bundle: D3 v7.9.0 + Highcharts (modern)
- ✅ Microscope pages: D3 v3 + NVD3 (legacy, via CDN)
- ✅ No conflicts: Legacy libs loaded as globals, main bundle imports modern versions
- ✅ No duplication: Only one `microscope.js` file (this one)

## Future Work

When time permits, `microscope.js` should be migrated to use the spectra
viewer in `packages/spectra`.

See git history for previous migration attempt and complexity analysis.

## Testing

End-to-end tests exist in `backend/fpbase/tests/test_end2end.py`:

- `test_microscopes()` - Tests microscope detail pages
- `test_embedscope()` - Tests embeddable microscope viewer

**Note**: These tests intentionally clear console logs instead of asserting no errors, as nvd3 can produce timing-related warnings that don't affect functionality.
