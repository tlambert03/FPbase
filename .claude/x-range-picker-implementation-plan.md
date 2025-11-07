# X-Axis Range Picker Implementation Plan

## Overview
Implement editable x-axis range inputs for the SpectraViewer that replace/overlay the min/max axis labels with interactive MUI TextField components.

## Key Requirements
1. Replace min/max x-axis labels with clickable MUI inputs
2. Bidirectional sync: inputs ↔ chart extremes ↔ zustand store
3. Individual min/max settable (user can set just one)
4. Persist to sessionStorage and URL params
5. Reset on "Remove All Spectra" button
6. Visual consistency with existing chart labels

## Type Changes

### Update ChartOptions type
File: `packages/spectra/src/types/index.ts`

```typescript
export interface ChartOptions {
  // ... existing fields
  extremes: [number | null, number | null] | null  // UPDATED
  // ... rest
}
```

**Rationale**:
- Outer `null` = no zoom applied
- Inner `null` values = individual extreme can be unset
- Examples:
  - `null` → No zoom (full range)
  - `[400, 700]` → Zoom from 400-700nm
  - `[null, 700]` → Zoom max only (min autoscales)
  - `[400, null]` → Zoom min only (max autoscales)

## Architecture

### Component Structure
```
XAxisRangeInputs (NEW component)
├─ MinInput (MUI TextField)
├─ MaxInput (MUI TextField)
└─ ClearButton (IconButton with CloseIcon)
```

**Key Design Decisions**:
- Component must be rendered INSIDE the XAxis component tree to use `useAxis()` hook
- Use absolute positioning to overlay inputs at bottom of chart
- Local state for input values (controlled inputs with validation on blur)
- Ref-based flag to prevent infinite loops in bidirectional sync

## Event Flow

### Chart Zoom → Store Update
```
User drags to zoom
  → Highcharts zoom
  → afterSetExtremes event fires
  → Check if extremes changed
  → Update store via updateChartOptions()
  → Inputs re-render with new values
```

### Input Edit → Chart Update
```
User types in input
  → Local state updates (controlled input)
  → User blurs input (or presses Enter)
  → Validate input (number, optional)
  → Update store via updateChartOptions()
  → useEffect watches extremes
  → Call axis.setExtremes()
```

### Preventing Infinite Loops
- In `afterSetExtremes`: Compare current store extremes with event values, only update if changed
- When calling `axis.setExtremes()`: Don't need flag because comparison in handler prevents loop

## Implementation Details

### 1. Update ChartOptions Type
File: `packages/spectra/src/types/index.ts`
- Change `extremes: [number, number] | null` to `extremes: [number | null, number | null] | null`

### 2. Update clearAllSpectra
File: `packages/spectra/src/store/spectraStore.ts`

```typescript
clearAllSpectra: () =>
  set((state) => ({
    activeSpectra: [],
    activeOverlaps: [],
    hiddenSpectra: [],
    exNorm: defaults.exNorm,
    chartOptions: {
      ...state.chartOptions,
      extremes: null  // ADD THIS
    },
    customFilters: {},
    customLasers: {},
    overlapCache: {},
  })),
```

### 3. Create XAxisRangeInputs Component
File: `packages/spectra/src/Components/SpectraViewer/XAxisRangeInputs.tsx` (NEW)

Key features:
- Use `useAxis()` to get axis instance
- Local state for input values (controlled)
- Validation on blur/Enter (not on change)
- Clear button when extremes are set
- Absolute positioning at bottom of chart

### 4. Update SpectraViewer
File: `packages/spectra/src/Components/SpectraViewer/SpectraViewer.jsx`

Changes needed:
- Update xAxis.events.afterSetExtremes to sync with store
- Integrate XAxisRangeInputs component (must be child of XAxis)
- Handle the special case where chartOptions.zoomType is set (line 144-149)

### 5. URL Params (NO CHANGES NEEDED)
File: `packages/spectra/src/utils/urlParams.ts`

Already handles individual extremes:
- Parse: Reads xMin/xMax separately
- Serialize: Writes each if not null
- Both support partial extremes

## Event Handler Pattern

### afterSetExtremes Handler
```typescript
afterSetExtremes: function(event) {
  const { min, max, userMin, userMax } = event

  // Reset case
  if (min === null && max === null) {
    if (chartOptions.extremes !== null) {
      updateChartOptions({ extremes: null })
    }
    return
  }

  // Zoom case - use userMin/userMax (actual zoom values)
  const newMin = userMin ?? min
  const newMax = userMax ?? max

  // Only update if changed (prevents loops)
  if (chartOptions.extremes?.[0] !== newMin ||
      chartOptions.extremes?.[1] !== newMax) {
    updateChartOptions({ extremes: [newMin ?? null, newMax ?? null] })
  }
}
```

## Visual Design

### Input Styling
- Size: `size="small"` (MUI)
- Width: 80px each
- Background: `rgba(255, 255, 255, 0.9)`
- Type: `type="number"` for numeric keyboard
- Font size: `0.75rem` (match axis labels)
- Padding: `4px 8px`
- Text align: center

### Positioning
```css
position: absolute;
bottom: 2px;
left: 10px;
right: 10px;
display: flex;
align-items: center;
gap: 8px;
z-index: 100;
pointer-events: none;  /* Allow clicks through to chart */
```

Individual inputs/button need `pointerEvents: 'auto'`

## Edge Cases

1. **Invalid input**: On blur, validate and revert to previous value if invalid
2. **Empty input**: Treat as null (unset that extreme)
3. **Reset zoom button**: Highcharts sends min=null, max=null
4. **Individual extremes**: Support via [null, 700] or [400, null]
5. **URL with only xMin or xMax**: Already supported by parsing logic
6. **Initialization**: Don't sync on mount if extremes already match

## Testing Checklist

- [ ] User zooms chart → inputs update
- [ ] User edits min input → chart zooms
- [ ] User edits max input → chart zooms
- [ ] User edits both → chart zooms to range
- [ ] User clears input → that extreme resets
- [ ] Chart reset button → inputs clear
- [ ] "Remove all spectra" button → inputs clear + chart resets
- [ ] URL with xMin/xMax → inputs populate
- [ ] URL with only xMin → min populates, max autoscales
- [ ] URL with only xMax → max populates, min autoscales
- [ ] Share button → URL includes current zoom
- [ ] Session storage → zoom persists on reload
- [ ] Invalid input → reverts to previous value
- [ ] Empty input → clears that extreme
- [ ] Enter key → commits input
- [ ] Chart zoom when inputs hidden (showX=false) → still works

## Files to Modify

1. ✅ `packages/spectra/src/types/index.ts` - Update ChartOptions.extremes type
2. ✅ `packages/spectra/src/Components/SpectraViewer/XAxisRangeInputs.tsx` - NEW component
3. ✅ `packages/spectra/src/Components/SpectraViewer/SpectraViewer.jsx` - Integrate inputs
4. ✅ `packages/spectra/src/store/spectraStore.ts` - Update clearAllSpectra

## Files NOT Modified (already correct)

- `packages/spectra/src/utils/urlParams.ts` - Already handles individual extremes
- `packages/spectra/src/defaults.ts` - extremes: null is correct default
- Session storage - Handled by Zustand persist middleware

## Implementation Order

1. Update types (ChartOptions)
2. Update clearAllSpectra to reset extremes
3. Create XAxisRangeInputs component
4. Integrate into SpectraViewer
5. Test all scenarios

## Notes

- react-jsx-highcharts requires components to be children of the axis to use useAxis()
- Need to handle the case where chartOptions.zoomType is set (line 144-149 in SpectraViewer)
- Comparison-based loop prevention is simpler than ref-based flags
- URL params already support individual extremes - no changes needed!
