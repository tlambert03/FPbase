# Migration Strategy: Full History Preservation (No Squashing)

## Overview

Move all 58 old migrations back, then add incremental migrations to transform the schema while maintaining ability to migrate any historical backup.

## File Structure Changes

**Before:**

```
migrations/
├── 0001_initial.py (NEW schema, complete)
├── 0002_update_ocfluoreff_to_use_direct_fk.py
migrations_old/
├── 0001_initial.py through 0058_*.py (OLD schema)
```

**After:**

```
migrations/
├── 0001_initial.py (moved from migrations_old)
├── 0002_* through 0058_* (moved from migrations_old)
├── 0059_add_fluorophore_and_new_models.py (NEW)
├── 0060_convert_state_to_mti.py (NEW)
├── 0061_migrate_dye_to_container_and_dyestate.py (NEW)
├── 0062_migrate_data.py (NEW - RunPython)
├── 0063_update_foreign_keys.py (NEW)
├── 0064_cleanup_old_fields.py (NEW)
```

## Migration Sequence

### Migration 0059: Add Fluorophore and New Models

**Purpose:** Create new schema elements non-destructively

**Operations:**

1. Create `Fluorophore` model (standalone for now)
2. Create `FluorescenceMeasurement` model with FK to Fluorophore
3. Create `DyeState` model (MTI child of Fluorophore)
4. Create NEW `Dye` container model with temp table name `proteins_dye_new`
5. Add `owner_fluor` FK to Spectrum (nullable, for transition)
6. Add `fluor` FK to OcFluorEff (nullable, for transition)

**Key:** All additions, no deletions. Old schema still works.

### Migration 0060: Convert State to MTI Child

**Purpose:** Make existing State model inherit from Fluorophore

**Operations:**

1. Add `fluorophore_ptr` OneToOneField to State (nullable initially)
2. For each existing State record:
   - Create corresponding Fluorophore record
   - Set State.fluorophore_ptr
3. Make fluorophore_ptr non-nullable
4. Add MTI meta options to State

**Challenge:** This is complex in Django. Alternative approach:

- Create StateNew with proper MTI
- Copy data State → StateNew  
- Rename tables later

**Decision needed:** Which approach for State MTI?

### Migration 0061: Migrate Dye → Container + DyeState

**Purpose:** Split old Dye into container Dye + DyeState

**RunPython Operations:**

```python
def migrate_dyes(apps, schema_editor):
    OldDye = apps.get_model('proteins', 'Dye')  # Uses old table
    NewDye = apps.get_model('proteins', 'DyeNew')  # Temp container
    DyeState = apps.get_model('proteins', 'DyeState')
    
    for old_dye in OldDye.objects.all():
        # Create container (copy all non-fluorescence fields)
        new_dye = NewDye.objects.create(
            name=old_dye.name,
            slug=old_dye.slug,
            # Copy: description, aliases, references, etc.
        )
        
        # Create DyeState (one per old Dye)
        # Fluorophore parent auto-created by MTI
        DyeState.objects.create(
            dye=new_dye,
            # Fluorophore fields:
            label=old_dye.name,
            slug=old_dye.slug,
            entity_type='dye',
        )
```

### Migration 0062: Migrate Fluorescence Data

**Purpose:** Create FluorescenceMeasurement records from old State/Dye data

**RunPython Operations:**

```python
def migrate_fluorescence_data(apps, schema_editor):
    State = apps.get_model('proteins', 'State')
    DyeState = apps.get_model('proteins', 'DyeState')
    FluorescenceMeasurement = apps.get_model('proteins', 'FluorescenceMeasurement')
    
    # Migrate State fluorescence data
    for state in State.objects.all():
        # Get reference (may be None)
        ref = state.protein.primary_reference if state.protein else None
        
        # Only create if has fluorescence data
        if state.ex_max or state.em_max or state.qy:
            FluorescenceMeasurement.objects.create(
                fluorophore=state.fluorophore_ptr,
                reference=ref,
                ex_max=state.ex_max,
                em_max=state.em_max,
                em_std=state.em_std,
                ext_coeff=state.ext_coeff,
                qy=state.qy,
                brightness=state.brightness,
                pka=state.pka,
                lifetime=state.lifetime,
                # ... all other fluorescence fields
            )
            # save() triggers rebuild_attributes() → materializes to Fluorophore
    
    # Migrate DyeState fluorescence data
    # (similar pattern)
```

### Migration 0063: Update Foreign Keys

**Purpose:** Point all FKs to new schema

**RunPython Operations:**

```python
def update_spectrum_ownership(apps, schema_editor):
    Spectrum = apps.get_model('proteins', 'Spectrum')
    State = apps.get_model('proteins', 'State')
    DyeState = apps.get_model('proteins', 'DyeState')
    
    # Update spectra owned by States
    for spectrum in Spectrum.objects.filter(owner_state__isnull=False):
        state = State.objects.get(id=spectrum.owner_state_id)
        spectrum.owner_fluor = state.fluorophore_ptr
        spectrum.save(update_fields=['owner_fluor'])
    
    # Update spectra owned by old Dyes
    # Need to find corresponding DyeState by slug
    for spectrum in Spectrum.objects.filter(owner_dye__isnull=False):
        old_dye = spectrum.owner_dye
        dyestate = DyeState.objects.get(slug=old_dye.slug)
        spectrum.owner_fluor = dyestate.fluorophore_ptr
        spectrum.save(update_fields=['owner_fluor'])

def update_ocfluoreff(apps, schema_editor):
    # Similar pattern - already in current 0002, adapt to 0063
```

### Migration 0064: Cleanup Old Schema

**Purpose:** Remove old fields and tables

**Operations:**

1. Drop fields from State:
   - All fluorescence property fields (ex_max, em_max, qy, etc.)
   - Keep: protein FK, name, maturation, etc.

2. Drop fields from Spectrum:
   - owner_state FK
   - owner_dye FK
   - Keep: owner_fluor FK (make non-nullable)

3. Rename Dye tables:
   - Drop old proteins_dye table
   - Rename proteins_dye_new → proteins_dye

4. Drop any other deprecated fields/tables

## Testing Protocol

### Test 1: Fresh Install (New Database)

```bash
dropdb fpbase && createdb fpbase
uv run python backend/manage.py migrate
uv run pytest --create-db
```

**Expected:** ✅ Runs 0001-0064, all tests pass

### Test 2: Production Migration Simulation

```bash
just pgpull  # Drops local, pulls production, runs migrate
```

**Expected:**

- Detects migrations 0001-0058 already applied
- Runs only 0059-0064
- ✅ Data migrated correctly

### Test 3: Partial Backup Migration

```bash
# Restore backup from middle of history (e.g., only 0001-0034 applied)
uv run python backend/manage.py migrate
```

**Expected:**

- Runs 0035-0064 in sequence
- ✅ Migrates successfully

### Test 4: Data Integrity Verification

```python
# In shell_plus after migration
# Verify State MTI
assert State.objects.filter(fluorophore_ptr__isnull=True).count() == 0

# Verify DyeState created
old_dye_count = 50  # Note before migration
assert DyeState.objects.count() == old_dye_count

# Verify Spectrum ownership migrated
assert Spectrum.objects.filter(owner_fluor__isnull=True).count() == 0
assert Spectrum.objects.filter(owner_state__isnull=False).count() == 0

# Verify Fluorophore materialized data
for state in State.objects.all():
    if state.protein.primary_reference:
        assert state.fluorophore_ptr.ex_max is not None
```

## Key Decision Points

**Q1: How to handle State MTI conversion?**

- Option A: Add fluorophore_ptr to existing State, populate, convert
- Option B: Create StateNew with MTI, migrate data, rename table
- **Recommendation:** Option A if possible (cleaner), Option B if Django doesn't support it

**Q2: Order of data migration vs FK updates?**

- Must create all Fluorophore/DyeState records BEFORE updating Spectrum FKs
- Order: 0059 (tables) → 0060 (State MTI) → 0061 (DyeState) → 0062 (data) → 0063 (FKs) → 0064 (cleanup)

**Q3: Handle BleachMeasurement, OSERMeasurement FKs?**

- These reference State via FK
- Need to verify they still work after State becomes MTI child
- May need FK updates in 0063

## Implementation Steps

1. **Move migrations:** `mv migrations_old/* migrations/`
2. **Delete current:** `rm migrations/0001_initial.py migrations/0002_*.py`
3. **Create 0059-0064** with schema operations above
4. **Test with fresh DB:** Verify all 64 migrations run
5. **Test with pgpull:** Verify production migration works
6. **Write verification script:** Data integrity checks
7. **Deploy to production**
8. **Delete migrations_old/** after success

## Risk Mitigation

- ✅ Full history preserved forever
- ✅ Any backup can migrate
- ✅ Incremental migrations easier to debug than big-bang
- ✅ Can test each migration step independently
- ⚠️ More migrations = more complexity
- ⚠️ MTI conversion is tricky, needs careful testing
