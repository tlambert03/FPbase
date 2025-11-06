# Icon System

FPbase uses an internal SVG icon library for optimal performance. Icons are rendered as inline `<svg>` elements, eliminating the need for Font Awesome's web fonts.

## Why This Approach?

- **Smaller bundle**: ~40 KB vs 200-800 KB with Font Awesome CDN
- **No FOUC**: Icons render immediately (no Flash of Unstyled Content)
- **SSR-friendly**: Works without JavaScript
- **Single source of truth**: Easy to swap icon libraries (e.g., migrate to Lucide)

## Django Templates

### Basic Usage

```django
{% load icon_tags %}

{# Simple icon #}
{% icon 'trash' %}

{# With CSS classes #}
{% icon 'edit' class='text-primary fa-lg' %}

{# With custom attributes #}
{% icon 'spinner' class='fa-spin' aria_label='Loading' %}
```

### Available Icons

All icons accept kebab-case names:
- `trash`, `edit`, `search`, `download`, `upload`
- `check`, `times`, `plus`, `minus`, `info-circle`
- `arrow-left`, `arrow-right`, `cog`, `user`, `heart`
- `github`, `google`, `orcid`, `x-twitter` (brands)

See full list: `backend/fpbase/icons.py` → `ICON_ALIASES`

### Size Modifiers

Use Font Awesome's size classes:

```django
{% icon 'home' class='fa-lg' %}     {# 1.33x larger #}
{% icon 'home' class='fa-2x' %}     {# 2x larger #}
{% icon 'home' class='fa-3x' %}     {# 3x larger #}
```

### Animations

```django
{% icon 'spinner' class='fa-spin' %}    {# Smooth rotation #}
{% icon 'cog' class='fa-pulse' %}       {# Stepped rotation #}
```

### In Buttons

```django
<button class="btn btn-danger">
  {% icon 'trash' %} Delete
</button>

<a href="#" class="btn btn-primary">
  {% icon 'plus' class='mr-2' %} Add New
</a>
```

## React/TypeScript

For React components, import icons directly:

```typescript
import { fasTrash, fasEdit } from '@/icons/fa-icons'

function Icon({ icon, className = '' }) {
  return (
    <svg
      viewBox={`0 0 ${icon.width} ${icon.height}`}
      fill="currentColor"
      className={className}
    >
      <path d={icon.path} />
    </svg>
  )
}

// Usage
<button><Icon icon={fasTrash} className="w-4 h-4" /> Delete</button>
```

## Adding New Icons

1. Update icon list in `scripts/extract-fa-icons.mjs`:
   ```javascript
   const iconMap = {
     // Add your icon here
     'my-new-icon': { pkg: 'solid', name: 'faMyNewIcon' },
   }
   ```

2. Regenerate icon files:
   ```bash
   node scripts/extract-fa-icons.mjs
   ```

3. Use in templates:
   ```django
   {% icon 'my-new-icon' %}
   ```

## Migrating from Font Awesome `<i>` Tags

### Old Way (Font Awesome CDN)
```html
<i class="fa fa-trash"></i>
<i class="fa fa-edit fa-lg"></i>
```

### New Way (Internal SVG)
```django
{% load icon_tags %}
{% icon 'trash' %}
{% icon 'edit' class='fa-lg' %}
```

## Switching to Lucide (Future)

When ready to migrate to Lucide, you only need to update:
1. `scripts/extract-fa-icons.mjs` → `scripts/extract-lucide-icons.mjs`
2. `backend/fpbase/templatetags/icon_tags.py` (no changes needed!)
3. `frontend/src/icons/*` files

Templates using `{% icon 'name' %}` will continue working!

## File Structure

```
backend/fpbase/
  icons.py                          # Generated: Icon data (Python)
  templatetags/
    icon_tags.py                    # Template tag implementation

frontend/src/icons/
  fa-icons.ts                       # Generated: Icon data (TypeScript)
  icon-replacer.js                  # Legacy: JS replacement (optional)
  icons.scss                        # Icon styles

scripts/
  extract-fa-icons.mjs              # Icon extraction script
```

## Performance

- **Before**: 200-800 KB (Font Awesome fonts + CSS)
- **After**: ~40 KB (inline SVG data, tree-shakeable)
- **Savings**: 4-20x smaller bundle size

## Troubleshooting

### Icon not found
```django
<!-- Icon "unknown-icon" not found -->
```
**Solution**: Check available icons in `backend/fpbase/icons.py` → `ICON_ALIASES`

### Icons not styled correctly
**Solution**: Ensure icon styles are loaded (`frontend/src/icons/icons.scss`)

### Need both kebab-case and camelCase?
Both work! `{% icon 'trash' %}` and `{% icon 'fasTrash' %}` are equivalent.
