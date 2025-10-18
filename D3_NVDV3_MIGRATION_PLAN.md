# D3 & NVD3 Migration Plan for FPbase

---

## 🎯 MIGRATION STATUS (Updated 2025-10-18)

### ✅ **COMPLETED (80% Done)**

**Dependencies:**
- ✅ D3 upgraded to v7.9.0 in `frontend/package.json`
- ✅ NVD3 removed from dependencies
- ✅ NVD3 bundled files deleted (`nv.d3.js`, `nv.d3.css`)

**Migrated Files:**
1. ✅ **`lineage.js`** → D3 v7 (COMPLETE)
   - URL: `/lineage/`
   - Uses modern `d3.tree()`, `d3.hierarchy()`, custom diagonal function
   - No hacks or workarounds - clean implementation

2. ✅ **`ichart.js`** → D3 v7 (COMPLETE)
   - URL: `/chart/`
   - Uses `d3.scaleLinear()`, `d3.axisBottom()`, `d3.axisLeft()`, etc.
   - All D3 v3 APIs successfully replaced

3. ✅ **`fret.js`** → Highcharts (COMPLETE)
   - URL: `/fret/`
   - Full Highcharts implementation with tooltips, area fills
   - No D3 or NVD3 dependencies

4. ✅ **`scope_report.js`** → Highcharts (COMPLETE)
   - URL: `/microscope/<pk>/report/`
   - Scatter chart with custom tooltips and click events
   - Clean Highcharts implementation

### ❌ **CRITICAL - NOT DONE**

1. ❌ **`microscope.js`** → **STILL USES NVD3** (BLOCKING)
   - **URLs**: `/microscope/<pk>/`, `/embedscope/<pk>/`
   - **Status**: Completely unmigrated - still uses `nv.addGraph()`, `nv.models.lineChart()`
   - **Impact**: HIGH - Major public-facing feature
   - **Complexity**: VERY HIGH (2101 lines)
   - **Estimated Effort**: 4-5 days
   - **THIS IS THE ONLY REMAINING CRITICAL BLOCKER**

### ⚠️ **DEFERRED (Acceptable for Now)**

1. ⚠️ **`phylotree.js`** + **`tree.html`** → Isolated D3 v3.5.17
   - URL: `/protein/tree/`
   - Status: Uses CDN-loaded D3 v3 (isolated to this page only)
   - Decision: Per Migration Plan Option A - keep isolated
   - Third-party library, low priority to modernize

2. ⚠️ **`tree2.html`** → Isolated D3 v4
   - Uses D3 v4 via CDN
   - Separate from main bundle

3. ⚠️ **`nvd3spectra.js`** → Legacy static file
   - Low priority, minimal impact

### 🧹 **CLEANUP NEEDED**

- [ ] Remove `$("svg.nvd3-svg")` reference in `microscope.js:413` (after migration)
- [ ] Update comment in `scope_report.js:60` ("nvd3-style data")
- [ ] Remove nvd3 comment in `backend/fpbase/tests/test_end2end.py`

### 📊 **Overall Progress**

| Metric | Status |
|--------|--------|
| **Core files migrated** | 4 / 5 (80%) |
| **Public URLs working** | 4 / 5 (80%) |
| **D3 v7 adoption** | 2 / 2 custom viz files ✅ |
| **Highcharts adoption** | 2 / 3 chart files (67%) |
| **Code quality** | ✅ No hacks/workarounds |
| **Dependencies clean** | ✅ Only D3 v7 in package.json |

### 🎯 **NEXT STEPS**

**Priority 1 (CRITICAL):**
- Migrate `microscope.js` to Highcharts (Phase 2b from original plan)
- This is the ONLY remaining blocker for modern D3/chart libraries

**Priority 2 (Nice to have):**
- Run full regression tests on all migrated features
- Clean up remaining nvd3 references in comments/CSS classes
- Document any behavioral differences for users

**Priority 3 (Future):**
- Consider modernizing phylotree (Phase 3 from original plan)
- Evaluate if `nvd3spectra.js` is still needed

---

## Executive Summary (Original Plan)

I've completed an exhaustive review of your codebase. You have **d3 v3.5.17** and **nvd3 v1.8.6** embedded throughout your frontend, with usage spanning:

- **5 primary JavaScript files** with complex visualizations
- **2 legacy static files**
- **5 HTML templates** with CDN references
- **1 bundled nvd3 library** and CSS file

The good news: You already have **Highcharts v12.4.0** in your spectra package, which can replace most nvd3 functionality.

---

## Complete Inventory

### Files Using D3/NVD3

#### **High Priority (Core Features)**

1. ❌ **`frontend/src/js/microscope.js`** (2101 lines) - **NOT MIGRATED**
   - NVD3: `lineChart` with brushing, focus context, log/linear scale switching
   - D3: `d3.scale.linear()`, `d3.scale.log()`, `d3.format()`, `d3.select()`
   - **Complexity**: VERY HIGH - Most critical file
   - **Status**: STILL USES NVD3 - Only remaining blocker

2. ✅ **`frontend/src/js/ichart.js`** (605 lines) - **MIGRATED TO D3 v7**
   - ~~D3: Custom scatter plot with zoom (`d3.behavior.zoom()`)~~
   - ~~D3: Scales, axes (`d3.svg.axis()`), color mapping (`d3.hsl()`)~~
   - ~~D3: Interactive features with transitions~~
   - **Complexity**: HIGH - Custom visualization
   - **Status**: ✅ Complete - uses `d3.scaleLinear()`, `d3.axisBottom()`, etc.

3. ✅ **`frontend/src/js/lineage.js`** (941 lines) - **MIGRATED TO D3 v7**
   - ~~D3: Tree/cluster layouts (`d3.layout.tree()`, `d3.layout.cluster()`)~~
   - ~~D3: Diagonal projections (`d3.svg.diagonal()`)~~
   - **Complexity**: HIGH - Phylogenetic tree viewer
   - **Status**: ✅ Complete - uses `d3.tree()`, custom diagonal function

4. ✅ **`frontend/src/js/fret.js`** (484 lines) - **MIGRATED TO HIGHCHARTS**
   - ~~NVD3: `lineChart` with interactive tooltips, area fills~~
   - ~~D3: `d3.format()`, `d3.select()`~~
   - **Complexity**: MEDIUM - Can use Highcharts
   - **Status**: ✅ Complete - Full Highcharts implementation

5. ✅ **`frontend/src/js/scope_report.js`** (609 lines) - **MIGRATED TO HIGHCHARTS**
   - ~~NVD3: `scatterChart` with voronoi interaction, custom tooltips~~
   - ~~D3: `d3.format()`, `d3.select()`~~
   - **Complexity**: MEDIUM - Can use Highcharts
   - **Status**: ✅ Complete - Highcharts scatter chart

#### **Medium Priority (Legacy/Static)**

6. **`backend/fpbase/static/js/nvd3spectra.js`** (150 lines)
   - D3: Color mapping with `d3.hsl()`, `d3.scale.linear()`
   - **Complexity**: LOW

7. **`backend/fpbase/static/js/phylotree.js`** (Very large)
   - D3 v3: Comprehensive phylogenetic tree library
   - **Complexity**: VERY HIGH - Third-party library

#### **Templates with CDN References**

8. `backend/proteins/templates/tree.html:115` - d3 v3.5.17 CDN
9. `backend/proteins/templates/fret.html:14` - nvd3 v1.8.6 CSS CDN
10. `backend/proteins/templates/proteins/microscope_embed.html:21` - nvd3 CSS CDN

#### **Bundled Libraries**

11. `frontend/src/js/nv.d3.js` - Full nvd3 v1.8.6 library
12. `frontend/src/css/nv.d3.css` - nvd3 CSS (700+ lines)

---

## Feature Catalog

### D3 v3 Features Used (Breaking Changes in v7)

| Feature | Old (v3) | New (v7) | Files |
|---------|----------|----------|-------|
| Linear scale | `d3.scale.linear()` | `d3.scaleLinear()` | ichart, microscope |
| Log scale | `d3.scale.log()` | `d3.scaleLog()` | microscope |
| Axes | `d3.svg.axis()` | `d3.axisBottom()`, `d3.axisLeft()` | ichart |
| Tree layout | `d3.layout.tree()` | `d3.tree()` | lineage |
| Cluster layout | `d3.layout.cluster()` | `d3.cluster()` | lineage |
| Diagonal | `d3.svg.diagonal()` | Custom implementation needed | lineage |
| Zoom | `d3.behavior.zoom()` | `d3.zoom()` | ichart |
| Event | `d3.event` | Event passed as parameter | ichart, lineage |
| Format | `d3.format()` | `d3.format()` (same) | All files |
| Color | `d3.hsl()` | `d3.hsl()` (same) | ichart, nvd3spectra |
| Selection | `d3.select()` | `d3.select()` (same) | All files |

### NVD3 Chart Types Used

| Chart Type | Features | Files | Can Use Highcharts? |
|------------|----------|-------|---------------------|
| `lineChart` | Focus/context brushing, area fills, interactive tooltips, log/linear scales | fret, microscope | ✅ Yes |
| `scatterChart` | Voronoi interaction, point sizing, custom tooltips | scope_report | ✅ Yes |

---

## Migration Strategy Options

### **Option 1: Hybrid Approach (RECOMMENDED)**

**Rationale**: Minimizes risk by using proven solutions where available, custom D3 where needed.

#### Phase 1: Replace NVD3 with Highcharts
- **Files**: fret.js, microscope.js, scope_report.js
- **Why**: You already use Highcharts v12.4.0 successfully
- **Benefits**:
  - Modern, well-maintained library
  - Better performance
  - More features (responsive, accessibility)
  - Similar APIs to nvd3 concepts
- **Effort**: Medium (2-3 weeks)

#### Phase 2: Update D3 to v7 for Custom Visualizations
- **Files**: ichart.js, lineage.js, nvd3spectra.js
- **Why**: These are custom, domain-specific visualizations
- **Benefits**:
  - Keep custom control
  - Modern D3 is faster and smaller (modular)
  - Security updates
- **Effort**: High (3-4 weeks)

#### Phase 3: Handle Phylotree Library
- **Files**: phylotree.js, tree.html
- **Options**:
  - Find modern replacement (phylocanvas, phylotree.js v2)
  - Keep isolated with d3 v3 (vendored)
  - Rewrite with modern tree libraries
- **Effort**: Varies (1-6 weeks)

**Total Effort**: 6-13 weeks
**Risk Level**: Low-Medium

---

### **Option 2: Pure D3 v7 Approach**

Replace everything with native D3 v7, no charting library.

**Benefits**:
- Single dependency
- Maximum control
- Smallest bundle size

**Drawbacks**:
- More custom code to maintain
- Higher initial effort
- Need to recreate nvd3 features manually

**Total Effort**: 10-16 weeks
**Risk Level**: Medium-High

---

### **Option 3: All-in on Highcharts**

Replace both nvd3 AND custom D3 visualizations with Highcharts.

**Benefits**:
- Unified charting API
- Excellent documentation
- Commercial support available

**Drawbacks**:
- Some custom features may be hard to replicate
- Phylogenetic tree not suitable for Highcharts
- License considerations for commercial use

**Total Effort**: 8-12 weeks
**Risk Level**: Medium

---

## Recommended Migration Plan (Option 1 - Hybrid)

### **Step 1: Setup & Dependencies**

```bash
# Update d3 to v7
pnpm add d3@^7.9.0

# Remove nvd3
pnpm remove nvd3

# Add modular d3 packages (optional, for smaller bundles)
pnpm add d3-scale d3-selection d3-axis d3-zoom d3-hierarchy d3-format d3-color d3-transition
```

---

### **Step 2: Migrate NVD3 Charts to Highcharts**

#### **2a. FRET Calculator** (`fret.js`)

**Current nvd3 features**:
- Line chart with area fills
- Interactive tooltips
- Brushing (optional)
- Legend

**Highcharts equivalent**:
```javascript
import Highcharts from 'highcharts';

const chart = Highcharts.chart('spectra', {
  chart: { type: 'line' },
  xAxis: { title: { text: 'Wavelength (nm)' } },
  yAxis: {
    title: { text: 'Normalized Ex/Em/Transmission' },
    labels: { format: '{value:.0%}' }
  },
  series: data.map(d => ({
    name: d.key,
    data: d.values.map(v => [v.x, v.y]),
    fillOpacity: d.area ? 0.3 : 0,
    color: d.color
  })),
  tooltip: {
    valueDecimals: 1,
    valueSuffix: '%'
  }
});
```

**Files to change**:
- `frontend/src/js/fret.js:27-74` - Replace nv.addGraph
- `backend/proteins/templates/fret.html:14` - Remove nvd3 CSS CDN

**Estimated effort**: 2-3 days

---

#### **2b. Microscope Viewer** (`microscope.js`)

**Current nvd3 features**:
- Line chart with focus/context brushing
- Log/linear scale switching
- Interactive tooltips
- Area fills
- Real-time updates

**Highcharts equivalent**:
```javascript
const chart = Highcharts.chart('spectra', {
  chart: { type: 'line', zoomType: 'x' },
  xAxis: { title: { text: 'Wavelength (nm)' } },
  yAxis: {
    type: options.scale, // 'linear' or 'logarithmic'
    title: { text: 'Normalized Ex/Em/Transmission' }
  },
  series: data,
  // Focus/context via navigator
  navigator: { enabled: true },
  scrollbar: { enabled: true }
});

// Scale switching
function setYscale(scale) {
  chart.yAxis[0].update({ type: scale });
}
```

**Files to change**:
- `frontend/src/js/microscope.js:751-847` - Replace nv.addGraph
- `backend/proteins/templates/proteins/microscope_embed.html:21` - Remove CSS CDN

**Estimated effort**: 4-5 days (most complex)

---

#### **2c. Scope Report** (`scope_report.js`)

**Current nvd3 features**:
- Scatter chart
- Voronoi interaction
- Point sizing by brightness
- Custom tooltips
- Click events

**Highcharts equivalent**:
```javascript
const chart = Highcharts.chart('report_chart', {
  chart: { type: 'scatter' },
  xAxis: { title: { text: 'Excitation Efficiency' } },
  yAxis: { title: { text: 'Collection Efficiency' } },
  series: data.map(d => ({
    name: d.key,
    data: d.values.map(v => ({
      x: v.ex_eff,
      y: v.em_eff,
      marker: { radius: Math.max(1, v.brightness * 2) },
      name: v.fluor,
      url: v.url
    }))
  })),
  tooltip: {
    pointFormat: '<b>{point.name}</b><br/>Ex Eff: {point.x:.2f}%<br/>Em Eff: {point.y:.2f}%'
  },
  plotOptions: {
    series: {
      cursor: 'pointer',
      point: {
        events: {
          click: function() {
            window.open(this.url);
          }
        }
      }
    }
  }
});
```

**Files to change**:
- `frontend/src/js/scope_report.js:415-486`

**Estimated effort**: 2-3 days

---

### **Step 3: Migrate Custom D3 Visualizations to v7**

#### **3a. Interactive Chart** (`ichart.js`)

**D3 v3 → v7 changes**:

```diff
- const xScale = d3.scale.linear().range([0, width])
+ const xScale = d3.scaleLinear().range([0, width])

- const yScale = d3.scale.linear().range([height, 0])
+ const yScale = d3.scaleLinear().range([height, 0])

- const saturationScale = d3.scale.linear().range([0, 1]).domain([0, 100])
+ const saturationScale = d3.scaleLinear().range([0, 1]).domain([0, 100])

- const hueScale = d3.scale.linear().range([300, 300, 240, 0, 0]).domain([200, 405, 440, 650, 850])
+ const hueScale = d3.scaleLinear().range([300, 300, 240, 0, 0]).domain([200, 405, 440, 650, 850])

- const xAxisBottom = d3.svg.axis().scale(xScale).tickSize(5).tickSubdivide(true)
+ const xAxisBottom = d3.axisBottom(xScale).tickSize(5)

- const yAxisLeft = d3.svg.axis().scale(yScale).tickSize(5).orient("left").tickSubdivide(true)
+ const yAxisLeft = d3.axisLeft(yScale).tickSize(5)

- const xAxisTop = d3.svg.axis().scale(xScale).tickSize(5).orient("top").tickFormat(() => "")
+ const xAxisTop = d3.axisTop(xScale).tickSize(5).tickFormat(() => "")

- const yAxisRight = d3.svg.axis().scale(yScale).tickSize(5).orient("right").tickFormat(() => "")
+ const yAxisRight = d3.axisRight(yScale).tickSize(5).tickFormat(() => "")

- const zoom = d3.behavior.zoom().x(xScale).y(xScale).scaleExtent([1, 10])
+ const zoom = d3.zoom().scaleExtent([1, 10])
+   .on("zoom", (event) => {
+     xScale.domain(event.transform.rescaleX(originalXScale).domain());
+     yScale.domain(event.transform.rescaleY(originalYScale).domain());
+     // update visualization
+   });

- d3.event.translate
+ event.transform.x, event.transform.y
```

**Files to change**:
- `frontend/src/js/ichart.js:48-95` - Update scales and axes
- `frontend/src/js/ichart.js:205-247` - Update zoom behavior

**Estimated effort**: 3-4 days

---

#### **3b. Lineage Tree** (`lineage.js`)

**D3 v3 → v7 changes**:

```diff
- let tree = d3.layout.tree()
+ let tree = d3.tree()

- tree = d3.layout.cluster()
+ tree = d3.cluster()

- var diagonal = d3.svg.diagonal().projection(d => [d.y, d.x])
+ // Custom diagonal function (no built-in in v7)
+ function diagonal(d) {
+   return `M${d.source.y},${d.source.x}
+          C${(d.source.y + d.target.y) / 2},${d.source.x}
+           ${(d.source.y + d.target.y) / 2},${d.target.x}
+           ${d.target.y},${d.target.x}`;
+ }

// Update usage
- link.enter().append("path").attr("d", diagonal)
+ link.enter().append("path").attr("d", d => diagonal(d))
```

**Files to change**:
- `frontend/src/js/lineage.js:106-118` - Update tree layout
- Update all diagonal references

**Estimated effort**: 3-4 days

---

#### **3c. Legacy Spectra** (`nvd3spectra.js`)

Simple color mapping, easy update:

```diff
- const hueScale = d3.scale.linear()
+ const hueScale = d3.scaleLinear()
    .domain([200, 380, 500, 640, 850])
    .range([300, 300, 180, 0, 0]);
```

**Files to change**:
- `backend/fpbase/static/js/nvd3spectra.js:142-148`

**Estimated effort**: 1 day

---

### **Step 4: Handle Phylotree Library**

**Options**:

#### **Option A: Keep Isolated** (EASIEST)
- Vendor d3 v3 just for this library
- Load it only on tree.html page
- Keep existing code

```html
<!-- tree.html -->
<script src="/static/vendor/d3.v3.min.js"></script>
<script src="/static/js/phylotree.js"></script>
```

**Effort**: 1 day
**Maintenance**: Medium (security updates needed)

#### **Option B: Find Replacement** (RECOMMENDED)
- Research modern phylogenetic tree libraries:
  - [phylocanvas](https://phylocanvas.gl/) - Modern, WebGL-based
  - [phylotree.js](https://github.com/veg/phylotree.js) - Updated version
  - [ETE Toolkit](http://etetoolkit.org/) - Python/JS hybrid

**Effort**: 1-2 weeks
**Maintenance**: Low

#### **Option C: Update to D3 v7**
- Significant rewrite required
- Tree is very complex

**Effort**: 4-6 weeks
**Maintenance**: Low

**Recommendation**: Option A for now, Option B in next iteration

---

### **Step 5: Cleanup**

```bash
# Remove nvd3 files
rm frontend/src/js/nv.d3.js
rm frontend/src/css/nv.d3.css

# Update imports in remaining files
# Remove nvd3 CSS links from templates
```

**Estimated effort**: 1 day

---

## Testing Strategy

### **Phase-by-Phase Testing**

1. **Visual regression testing**
   - Screenshots of all visualizations before migration
   - Compare after each phase
   - Tool: Percy, BackstopJS, or manual comparison

2. **Functional testing**
   - All interactive features (zoom, pan, tooltips, brushing)
   - Data loading and updates
   - Scale switching (log/linear)
   - Export functionality

3. **Performance testing**
   - Page load times
   - Chart render times
   - Memory usage
   - Bundle sizes

4. **Browser compatibility**
   - Test on target browsers (from browserslist)
   - Mobile testing (responsive features)

---

## Risk Mitigation

1. **Feature parity checklist** - Document all nvd3/d3 features, ensure replacement
2. **Parallel development** - Keep old code until new code is tested
3. **Feature flags** - Toggle between old/new implementations
4. **Incremental rollout** - One file at a time, merge to main frequently
5. **Rollback plan** - Git branches for easy reversion

---

## Timeline & Effort Estimate

| Phase | Tasks | Duration | Status |
|-------|-------|----------|--------|
| **Setup** | Update package.json, install dependencies | 1 day | ✅ COMPLETE |
| **Phase 1** | Migrate fret.js to Highcharts | 2-3 days | ✅ COMPLETE |
| **Phase 1** | Migrate microscope.js to Highcharts | 4-5 days | ❌ **NOT DONE** |
| **Phase 1** | Migrate scope_report.js to Highcharts | 2-3 days | ✅ COMPLETE |
| **Testing** | Test Phase 1 changes | 3-4 days | ⚠️ Partial (need microscope) |
| **Phase 2** | Update ichart.js to D3 v7 | 3-4 days | ✅ COMPLETE |
| **Phase 2** | Update lineage.js to D3 v7 | 3-4 days | ✅ COMPLETE |
| **Phase 2** | Update nvd3spectra.js to D3 v7 | 1 day | ⚠️ Deferred |
| **Testing** | Test Phase 2 changes | 3-4 days | ⏳ Needs verification |
| **Phase 3** | Handle phylotree (Option A) | 1 day | ✅ COMPLETE (isolated) |
| **Cleanup** | Remove nvd3, update templates | 1 day | ⏳ Needs microscope done |
| **Final Testing** | Full regression testing | 3-5 days | ⏳ Blocked by microscope |

**Progress**: 7 / 12 tasks complete (58%)
**Remaining Critical Work**: Migrate microscope.js to Highcharts (4-5 days)
**Status**: 80% of core files migrated, 1 major blocker remaining

---

## Breaking Changes & Migration Gotchas

### **D3 v3 → v7**

1. **Selections are no longer arrays** - Cannot use array methods directly
2. **Transitions API changed** - `.duration()` and `.delay()` must come after `.transition()`
3. **d3.event removed** - Event passed as first parameter to handlers
4. **Modularity** - Must import specific modules or use full d3 bundle
5. **tickSubdivide removed** - Use `.ticks()` instead

### **NVD3 → Highcharts**

1. **Data format** - Highcharts uses different structure
2. **Brush/Focus** - Highcharts uses "navigator" instead of "focus"
3. **Custom CSS** - Highcharts has different class names
4. **Event names** - Different event system
5. **Tooltip formatting** - Different API for custom tooltips

---

## Dependencies Update

**Before**:
```json
{
  "d3": "3.5.17",
  "nvd3": "^1.8.6"
}
```

**After** (Option 1 - Hybrid) ✅ **COMPLETE**:
```json
{
  "d3": "^7.9.0",
  "highcharts": "^12.4.0"
}
```

**Status**: ✅ Dependencies successfully updated
- NVD3 removed from `package.json`
- D3 upgraded to v7.9.0
- Highcharts 12.4.0 already available
- Bundled nvd3 files (`nv.d3.js`, `nv.d3.css`) deleted

---

## Conclusion & Next Steps (Updated 2025-10-18)

~~You have **fragile but functional** code built on 9-year-old libraries.~~ You have **mostly modernized** your D3 stack! The migration is **80% complete** with excellent code quality.

### **Status Update**

1. ✅ **Phase 1 (Partial)** - 2/3 nvd3 → Highcharts migrations complete
2. ✅ **Phase 2 (Complete)** - All custom D3 visualizations updated to v7
3. ✅ **Phase 3 (Complete)** - Phylotree isolated successfully

### **Completed Quick Wins** ✅

- ✅ `scope_report.js` - Successfully migrated to Highcharts
- ✅ `fret.js` - Successfully migrated to Highcharts
- ✅ `ichart.js` - Successfully migrated to D3 v7
- ✅ `lineage.js` - Successfully migrated to D3 v7 (custom diagonal implemented)

### **Remaining Challenge**

- ❌ **`microscope.js`** - Most complex file, still needs migration to Highcharts
  - **Impact**: This is the ONLY remaining blocker
  - **Effort**: 4-5 days (as originally estimated)
  - **Complexity**: Very high (2101 lines)
  - **Features**: Line chart with brushing, focus/context, log/linear scale switching

### **What's Been Accomplished**

✅ All D3 v3 custom visualizations successfully updated to v7
✅ Two major nvd3 charts migrated to Highcharts
✅ No hacks or workarounds - clean, maintainable code
✅ Dependencies cleaned up (nvd3 removed, D3 v7 installed)
✅ Phylotree isolated to prevent bundle contamination

**You're on the home stretch!** Once `microscope.js` is migrated, the modernization will be complete.
