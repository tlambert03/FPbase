🧪 Manual Testing Checklist

CRITICAL: Core Functionality (All Pages)

Test on any page first:

- [ ] Page loads without JavaScript errors (check browser console)
- [ ] Search autocomplete in navbar works (start typing a protein name)
- [ ] Bootstrap modals open/close properly
- [ ] Bootstrap tooltips/popovers display correctly
- [ ] No visual layout issues or missing styles


---

5. Protein Detail Page

URL: /protein/{slug}/ (e.g., /protein/egfp/)

Test:

- [ ] Page loads successfully
- [ ] If protein has PDB structure, LiteMol viewer loads
- [ ] SnapGene external resources section appears
- [ ] Lineage chart (if present) loads and renders
- [ ] Lineage chart scrolls to highlight current protein
- [ ] Spectra viewer (if present) displays correctly
- [ ] All interactive elements work

Changed: Multiple jQuery patterns → vanilla JS (scroll, lineage, events)

---

6. Protein Search

URL: /proteins/

Test:

- [ ] Search page loads
- [ ] Search functionality initializes
- [ ] Filter dropdowns work
- [ ] Search autocomplete works
- [ ] Results display correctly

Changed: $(function) wrapper → DOMContentLoaded

---

7. Protein Form (Add/Edit)

URL: /protein/create/ or /protein/{slug}/update/

Test:

- [ ] Form loads successfully
- [ ] "Add State" button appears and works (creates new formset)
- [ ] Delete button (×) on dynamic states works
- [ ] If updating existing protein, empty state is auto-removed on load
- [ ] When states count > 1, transition note appears
- [ ] When states count = 1, transition note hides
- [ ] Form validation works (try submitting with invalid name)
- [ ] Invalid form scrolls to top smoothly
- [ ] Form submission works

Changed: Extensive jQuery → vanilla JS conversion (formsets, validation, show/hide)

---

8. Bleach Measurement Forms

URL: /protein/{slug}/bleach/add/ or /bleach-comparison/{ref_id}/

Test:

- [ ] Form loads
- [ ] "Add Bleach Measurement" / "Add protein" button works
- [ ] Dynamic form rows can be added
- [ ] Delete button (×) on dynamic forms works
- [ ] If updating, empty form is auto-removed on load
- [ ] Form submission works

Changed: $(function) wrapper, dynamic form deletion → vanilla JS

---

9. Organism List

URL: /organisms/

Test:

- [ ] Page loads
- [ ] DataTable initializes (pagination, sorting, search)
- [ ] Table sorting works (click column headers)
- [ ] Table filtering/search works
- [ ] "Show 10/25/50/100/All" entries dropdown works
- [ ] Horizontal scroll works on mobile/narrow screens

Changed: $(function) → DOMContentLoaded for DataTables init

---

10. Organism Detail

URL: /organism/{id}/ (e.g., /organism/1/)

Test:

- [ ] Page loads
- [ ] Lineage chart (if present) loads via AJAX
- [ ] "Lineages derived from X" heading appears dynamically
- [ ] Lineage chart renders correctly
- [ ] Chart interactions work

Changed: jQuery DOM creation → vanilla JS

---

11. Lineage Page

URL: /lineage/{slug}/

Test:

- [ ] Page loads with loading spinner
- [ ] Loading message appears initially
- [ ] Lineage chart loads via AJAX
- [ ] Loading spinner disappears after load
- [ ] Lineage SVG renders correctly
- [ ] Search/toolbar features work
- [ ] Top scroll works

Changed: .each(), .html(), .empty() → vanilla JS

---

12. Old Spectra Viewer

URL: /spectra/ (if this route exists)

Test:

- [ ] Spectra visualization loads
- [ ] "Download Spectra" button works
- [ ] Clicking download triggers CSV export
- [ ] Only enabled (non-disabled) spectra are included

Changed: $("#downloadSpectra").click() → vanilla JS event listener

---

13. Spectrum Form

URL: /spectrum/create/ or /spectrum/{id}/update/

Test:

- [ ] Form loads with tabbed interface (File Upload / Manual Entry)
- [ ] Tabs switch correctly between file upload and manual data entry
- [ ] Category dropdown changes subtype options dynamically
- [ ] Owner type selection shows/hides protein vs non-protein fields
- [ ] File upload preview works
- [ ] Manual data entry preview works
- [ ] Preview shows correct chart/statistics
- [ ] "Edit Data" button returns to form
- [ ] Form validation prevents submission without data
- [ ] Submit button updates based on form state
- [ ] Form submission works

Changed: 3x $(document).ready() → DOMContentLoaded (complex form logic)

---

14. Microscope Report

URL: /microscope/{id}/

Test:

- [ ] Page loads
- [ ] ScopeReport widget initializes
- [ ] Select2 dropdown for probe filter works
- [ ] Probe filter is searchable
- [ ] DataTables loads correctly
- [ ] Report updates with selected filters

Changed: 2x jQuery wrappers → DOMContentLoaded

---

15. Microscope Embed

URL: /microscope/{id}/embed/

Test:

- [ ] Embedded microscope page loads
- [ ] Efficiency link has target="_blank" attribute
- [ ] External link icon appears next to link
- [ ] Link opens in new tab when clicked
- [ ] Embedded scope visualization works

Changed: .attr(), .append() → vanilla JS

---

16. Protein Collection Detail

URL: /collection/{slug}/

Test:

- [ ] Page loads
- [ ] Display type radio buttons appear (if multiple display options)
- [ ] Clicking display buttons switches between views
- [ ] Only selected display type is visible
- [ ] Other display types are hidden

Changed: Event handlers, show/hide → vanilla JS

---

17. Activity Page

URL: /activity/ (if exists, analytics page)

Test:

- [ ] Page loads
- [ ] Range selector dropdown appears
- [ ] Changing range selector shows/hides corresponding stats sections
- [ ] Correct stats section displays for selected range

Changed: Change handler, show/hide → vanilla JS

---

18. Reference Detail

URL: /reference/{id}/

Test:

- [ ] Page loads successfully
- [ ] "Add an excerpt" button appears
- [ ] Clicking button opens modal
- [ ] Modal appears ON TOP of backdrop (not behind)
- [ ] Modal can be closed via X button
- [ ] Modal can be closed by clicking outside
- [ ] Modal form submission works (if logged in)
- [ ] Info popover on "Excerpts" heading works

Changed: Removed z-index manipulation code (potential modal issue!)

---

19. Reference List

URL: /references/

Test:

- [ ] Page loads
- [ ] Any functionality on this page works

Changed: Minor change (check git diff if needed)

---

20. History/Revision Page

URL: /protein/{slug}/history/

Test:

- [ ] Page loads (staff only)
- [ ] Revision history displays
- [ ] "Revert to this revision" link works
- [ ] Clicking revert link submits form
- [ ] Page doesn't navigate away on link click (preventDefault works)

Changed: $(this).closest('form').submit() → vanilla JS

---
🎯 Priority Testing Order

1. Start here: Test one protein detail page (/protein/egfp/)
2. Then: Reference detail with modal (/reference/1/)
3. Then: Complex forms (protein form, spectrum form)
4. Then: Interactive features (chart, lineage, FRET)
5. Finally: Less critical pages (organism list, activity)

---
🚨 Watch For These Issues

- [ ] Console errors (open DevTools F12)
- [ ] Elements not appearing (things that should load dynamically)
- [ ] Modals appearing behind backdrop (especially reference excerpts)
- [ ] Forms not submitting or validating incorrectly
- [ ] Buttons not responding to clicks
- [ ] Animations missing or happening instantly
- [ ] Select2/DataTables not initializing
