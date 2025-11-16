# Reusable Button Components

This directory contains reusable Django template includes for Bootstrap 4 buttons used throughout FPbase.

## Files

- `_button.html` - Full-featured button component for all button types
- `_btn_circle.html` - Simplified component for circular icon-only buttons

## `_button.html` - Standard Button Component

### Basic Usage

```django
{% include 'includes/_button.html' with text="Submit" variant="primary" type="submit" %}
```

### Parameters

#### Required
- `text` - Button text content (can be empty for icon-only buttons)

#### Optional - Button Style
- `type` - Button type attribute (default: `"button"`)
  - Options: `"button"`, `"submit"`, `"reset"`
- `variant` - Bootstrap button variant (default: `"primary"`)
  - Options: `primary`, `secondary`, `info`, `danger`, `warning`, `success`, `light`, `dark`, `link`
- `outline` - Use outline variant (default: `False`)
- `size` - Button size (default: `""`)
  - Options: `sm`, `lg`, `""` (empty for default)
- `block` - Make button full width (default: `False`)
- `block_xs` - Apply `btn-block-xs` class (default: `False`)
- `circle` - Apply `btn-circle` class for icon-only circular buttons (default: `False`)

#### Optional - Icon
- `icon` - Icon name to include (default: `""`)
- `icon_class` - Additional classes for the icon (default: `""`)
- `icon_after` - Show icon after text instead of before (default: `False`)

#### Optional - Responsive Display
- `hide_sm_up` - Hide on sm screens and up, show only on xs (default: `False`)
  - Adds: `d-block d-sm-none`
- `hide_xs` - Hide on xs screens, show on sm and up (default: `False`)
  - Adds: `d-none d-sm-inline` (or `d-none d-sm-block` if `hide_xs_inline=False`)
- `hide_xs_inline` - When `hide_xs=True`, use inline display (default: `True`)

#### Optional - HTML Attributes
- `class` - Additional CSS classes (default: `""`)
- `id` - Button ID attribute (default: `""`)
- `disabled` - Disable the button (default: `False`)
- `attrs` - Additional HTML attributes as a string (default: `""`)

#### Optional - Data Attributes
- `data_toggle` - `data-toggle` attribute (default: `""`)
- `data_target` - `data-target` attribute (default: `""`)
- `data_dismiss` - `data-dismiss` attribute (default: `""`)

For other `data-*` attributes, use the `attrs` parameter.

### Examples

#### Primary submit button
```django
{% include 'includes/_button.html' with text="Submit" variant="primary" type="submit" %}
```
→ `<button type="submit" class="btn btn-primary">Submit</button>`

#### Small danger button
```django
{% include 'includes/_button.html' with text="Delete" variant="danger" size="sm" %}
```
→ `<button type="button" class="btn btn-danger btn-sm">Delete</button>`

#### Button with icon
```django
{% include 'includes/_button.html' with icon="add" text="Add Item" variant="info" %}
```
→ `<button type="button" class="btn btn-info">{% icon "add" %} Add Item</button>`

#### Icon after text
```django
{% include 'includes/_button.html' with text="External Link" icon="external-link" icon_after=True %}
```
→ `<button type="button" class="btn btn-primary">External Link {% icon "external-link" %}</button>`

#### Modal trigger button
```django
{% include 'includes/_button.html' with text="Open Modal" variant="info" data_toggle="modal" data_target="#myModal" %}
```
→ `<button type="button" class="btn btn-info" data-toggle="modal" data-target="#myModal">Open Modal</button>`

#### Outline button
```django
{% include 'includes/_button.html' with text="Cancel" variant="secondary" outline=True %}
```
→ `<button type="button" class="btn btn-outline-secondary">Cancel</button>`

#### Full-width block button
```django
{% include 'includes/_button.html' with text="Save" variant="primary" block=True %}
```
→ `<button type="button" class="btn btn-primary btn-block">Save</button>`

#### Mobile-only button (hidden on sm+)
```django
{% include 'includes/_button.html' with icon="sun" text="Edit States/Attributes" variant="info" block=True hide_sm_up=True class="mt-3" %}
```
→ `<button type="button" class="btn btn-info btn-block d-block d-sm-none mt-3">{% icon "sun" %} Edit States/Attributes</button>`

#### Desktop-only inline element (hidden on xs)
```django
{% include 'includes/_button.html' with text="Details" variant="link" size="sm" hide_xs=True %}
```
→ `<button type="button" class="btn btn-link btn-sm d-none d-sm-inline">Details</button>`

#### Button with custom attributes
```django
{% include 'includes/_button.html' with text="Custom" variant="info" attrs='data-custom="value" onclick="myFunction()"' %}
```
→ `<button type="button" class="btn btn-info" data-custom="value" onclick="myFunction()">Custom</button>`

#### Close button for modal
```django
{% include 'includes/_button.html' with text="Close" variant="secondary" data_dismiss="modal" %}
```
→ `<button type="button" class="btn btn-secondary" data-dismiss="modal">Close</button>`

## `_btn_circle.html` - Circular Icon Button Component

Simplified component specifically for icon-only circular buttons (common in tables and toolbars).

### Basic Usage

```django
{% include 'includes/_btn_circle.html' with icon="add" %}
```

### Parameters

#### Required
- `icon` - Icon name to display

#### Optional
- `variant` - Bootstrap button variant (default: `"info"`)
- `outline` - Use outline variant (default: `True`)
- `size` - Button size (default: `"sm"`)
- `type` - Button type attribute (default: `"button"`)
- `class` - Additional CSS classes
- `id` - Button ID attribute
- `disabled` - Disable the button
- `attrs` - Additional HTML attributes as a string
- `title` - Tooltip/title attribute
- `aria_label` - aria-label attribute
- `data_toggle`, `data_target`, `data_dismiss` - Common data attributes

### Examples

#### Basic circle button
```django
{% include 'includes/_btn_circle.html' with icon="add" %}
```
→ `<button type="button" class="btn btn-outline-info btn-sm btn-circle">{% icon "add" %}</button>`

#### Modal trigger circle button
```django
{% include 'includes/_btn_circle.html' with icon="info-i" data_toggle="modal" data_target="#fretInfoModal" %}
```
→ `<button type="button" class="btn btn-outline-info btn-sm btn-circle" data-toggle="modal" data-target="#fretInfoModal">{% icon "info-i" %}</button>`

#### Danger circle button with title
```django
{% include 'includes/_btn_circle.html' with icon="close" variant="danger" title="Remove item" %}
```
→ `<button type="button" class="btn btn-outline-danger btn-sm btn-circle" title="Remove item">{% icon "close" %}</button>`

#### Solid (non-outline) circle button
```django
{% include 'includes/_btn_circle.html' with icon="add" outline=False %}
```
→ `<button type="button" class="btn btn-info btn-sm btn-circle">{% icon "add" %}</button>`

## Migration Examples

Here are examples of migrating existing button code to use the new includes:

### Example 1: Submit button
**Before:**
```django
<button type="submit" class="btn btn-primary">Submit</button>
```

**After:**
```django
{% include 'includes/_button.html' with text="Submit" variant="primary" type="submit" %}
```

### Example 2: Button with icon and responsive classes
**Before:**
```django
<button type="button" class="btn btn-info btn-block mt-3 d-block d-sm-none">
    {% icon "sun" class_="mr-2" %}Edit States&sol;Attributes
</button>
```

**After:**
```django
{% include 'includes/_button.html' with icon="sun" text="Edit States&sol;Attributes" variant="info" block=True hide_sm_up=True class="mt-3" %}
```

### Example 3: Modal close button
**Before:**
```django
<button type="button" class="btn btn-secondary" data-dismiss="modal">Close</button>
```

**After:**
```django
{% include 'includes/_button.html' with text="Close" variant="secondary" data_dismiss="modal" %}
```

### Example 4: Small info button with modal trigger
**Before:**
```django
<button type="button" class='btn btn-sm btn-secondary mt-1' data-toggle="modal" data-target="#referenceModal">
    {% icon "collection" class_="mr-2" %}Add a reference
</button>
```

**After:**
```django
{% include 'includes/_button.html' with icon="collection" text="Add a reference" variant="secondary" size="sm" class="mt-1" data_toggle="modal" data_target="#referenceModal" %}
```

### Example 5: Circle button for comparison
**Before:**
```django
<button class="btn btn-sm btn-outline-info btn-circle comparison-btn"
        data-flash='1'
        data-action-url="{% url 'proteins:update-comparison' %}"
        data-object='{{protein.slug}}'
        data-op='add'>
    {% icon "add" %}
</button>
```

**After:**
```django
{% include 'includes/_btn_circle.html' with icon="add" class="comparison-btn" attrs='data-flash="1" data-action-url="{% url "proteins:update-comparison" %}" data-object="{{protein.slug}}" data-op="add"' %}
```

### Example 6: Info circle button with modal
**Before:**
```django
<button type="button" data-toggle="modal" data-target="#fretInfoModal" class='btn btn-info btn-circle ml-1'>
    {% icon "info-i" %}
</button>
```

**After:**
```django
{% include 'includes/_btn_circle.html' with icon="info-i" variant="info" outline=False class="ml-1" data_toggle="modal" data_target="#fretInfoModal" %}
```

### Example 7: Fetch button with input group
**Before:**
```django
<button type="submit" class="btn btn-info">Fetch</button>
```

**After:**
```django
{% include 'includes/_button.html' with text="Fetch" variant="info" type="submit" %}
```

### Example 8: Alert close button
**Before:**
```django
<button type="button" class="close" data-dismiss="alert" aria-label="Close">
    <span aria-hidden="true">&times;</span>
</button>
```

**After:**
This is a special Bootstrap close button that should remain as-is, or you could create a separate `_close_button.html` include if needed.

## Migration Strategy

1. **Start with high-frequency patterns**: Focus on the most common button patterns first (submit buttons, modal triggers, circle buttons)
2. **Test responsive behavior**: Pay special attention to buttons with `d-block d-sm-none` or `d-none d-sm-inline` classes
3. **Validate icon spacing**: The include automatically adds spacing between icons and text, so remove manual `mr-2` from icon classes
4. **Check data attributes**: Complex data attributes should be passed via the `attrs` parameter
5. **Test each migration**: After converting a button, verify the rendered HTML matches the original

## Notes

- The include automatically handles icon-text spacing, so you don't need `class_="mr-2"` on icons
- For buttons inside `<a>` tags, you may want to refactor to use the button include with appropriate onclick handlers or convert the anchor to a button
- The `btn-block-xs` class is custom and should be defined in your CSS if not already present
- Close buttons (`<button class="close">`) are a special Bootstrap pattern and may warrant a separate include if used frequently
