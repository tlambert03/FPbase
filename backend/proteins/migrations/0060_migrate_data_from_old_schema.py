# Generated manually for schema overhaul

import logging

from django.db import migrations

logger = logging.getLogger(__name__)


def migrate_state_data(apps, schema_editor):
    """Migrate State data from old schema to new Fluorophore + State MTI structure."""
    # Get models from migration state
    Fluorophore = apps.get_model("proteins", "Fluorophore")
    State = apps.get_model("proteins", "State")
    Protein = apps.get_model("proteins", "Protein")
    FluorescenceMeasurement = apps.get_model("proteins", "FluorescenceMeasurement")

    # Access old State data directly via raw SQL
    with schema_editor.connection.cursor() as cursor:
        cursor.execute("""
            SELECT id, created, modified, name, slug, is_dark,
                   ex_max, em_max, ext_coeff, qy, brightness,
                   lifetime, pka, twop_ex_max, "twop_peakGM", twop_qy,
                   maturation, protein_id
            FROM proteins_state_old
        """)

        for row in cursor.fetchall():
            (old_id, created, modified, name, slug, is_dark,
             ex_max, em_max, ext_coeff, qy, brightness,
             lifetime, pka, twop_ex_max, twop_peakgm, twop_qy,
             maturation, protein_id) = row

            # Get protein for label
            try:
                protein = Protein.objects.get(id=protein_id)
                label = f"{protein.name} ({name})" if name and name != "default" else protein.name
                # Handle empty/null slugs with guaranteed non-empty fallback
                if not slug or (isinstance(slug, str) and slug.strip() == ""):
                    # Try protein slug + name, or protein slug, or fallback to state ID
                    if name and name != "default":
                        state_slug = f"{protein.slug}-{name}" if protein.slug else f"state-{old_id}"
                    else:
                        state_slug = protein.slug if protein.slug else f"state-{old_id}"
                else:
                    state_slug = slug

                # Final safety check - ensure slug is not empty
                if not state_slug or state_slug.strip() == "":
                    state_slug = f"state-{old_id}"

            except Protein.DoesNotExist:
                print(f"Warning: Protein {protein_id} not found for State {old_id}, skipping")
                continue

            # Ensure slug uniqueness by checking if it already exists
            original_slug = state_slug
            base_slug = state_slug
            counter = 1
            while Fluorophore.objects.filter(slug=state_slug).exists():
                state_slug = f"{base_slug}-{counter}"
                counter += 1
                if counter > 100:
                    raise ValueError(
                        f"Could not generate unique slug for State {old_id} "
                        f"after 100 attempts (original: {original_slug})"
                    )

            if state_slug != original_slug:
                logger.warning(
                    f"Slug collision during State migration: {original_slug} -> {state_slug} "
                    f"(State ID: {old_id}, Protein: {protein.name})"
                )

            # Create Fluorophore parent (MTI will link automatically)
            fluorophore = Fluorophore.objects.create(
                created=created,
                modified=modified,
                label=label,
                slug=state_slug,
                entity_type="protein",
                ex_max=ex_max,
                em_max=em_max,
                ext_coeff=ext_coeff,
                qy=qy,
                brightness=brightness,
                lifetime=lifetime,
                pka=pka,
                twop_ex_max=twop_ex_max,
                twop_peakGM=twop_peakgm,  # Map from SQL result to model field
                twop_qy=twop_qy,
                is_dark=is_dark,
                
                
            )

            # Create State (MTI child) pointing to the Fluorophore
            # Use raw SQL to avoid Django MTI trying to update parent with empty values
            cursor.execute("""
                INSERT INTO proteins_state (fluorophore_ptr_id, name, protein_id, maturation)
                VALUES (%s, %s, %s, %s)
            """, [fluorophore.pk, name, protein_id, maturation])

            # Create FluorescenceMeasurement from old State data if there's any fluorescence data
            if any([ex_max, em_max, qy, ext_coeff, lifetime, pka]):
                # Get protein's primary reference (may be None)
                reference = protein.primary_reference if hasattr(protein, 'primary_reference') else None
                reference_id = protein.primary_reference_id if hasattr(protein, 'primary_reference_id') else None

                FluorescenceMeasurement.objects.create(
                    fluorophore=fluorophore,
                    reference_id=reference_id,
                    ex_max=ex_max,
                    em_max=em_max,
                    ext_coeff=ext_coeff,
                    qy=qy,
                    brightness=brightness,
                    lifetime=lifetime,
                    pka=pka,
                    twop_ex_max=twop_ex_max,
                    twop_peakGM=twop_peakgm,
                    twop_qy=twop_qy,
                    is_dark=is_dark,
                    is_trusted=True,  # Mark as trusted since it's the original data
                    
                    
                )

    print(f"Migrated {State.objects.count()} State records")


def migrate_dye_data(apps, schema_editor):
    """Migrate Dye data from old schema to new Dye container + DyeState structure."""
    Fluorophore = apps.get_model("proteins", "Fluorophore")
    Dye = apps.get_model("proteins", "Dye")
    DyeState = apps.get_model("proteins", "DyeState")
    FluorescenceMeasurement = apps.get_model("proteins", "FluorescenceMeasurement")

    # Access old Dye data via raw SQL
    with schema_editor.connection.cursor() as cursor:
        cursor.execute("""
            SELECT id, created, modified, name, slug, is_dark,
                   ex_max, em_max, ext_coeff, qy, brightness,
                   lifetime, pka, twop_ex_max, "twop_peakGM", twop_qy
            FROM proteins_dye_old
        """)

        for row in cursor.fetchall():
            (old_id, created, modified, name, slug, is_dark,
             ex_max, em_max, ext_coeff, qy, brightness,
             lifetime, pka, twop_ex_max, twop_peakgm, twop_qy) = row

            # Handle empty/null slugs
            if not slug or slug.strip() == "":
                dye_slug = f"dye-{old_id}"  # Use old ID as fallback
            else:
                dye_slug = slug

            # Ensure Dye slug uniqueness
            original_dye_slug = dye_slug
            base_dye_slug = dye_slug
            counter = 1
            while Dye.objects.filter(slug=dye_slug).exists():
                dye_slug = f"{base_dye_slug}-{counter}"
                counter += 1
                if counter > 100:
                    raise ValueError(
                        f"Could not generate unique slug for Dye {old_id} "
                        f"after 100 attempts (original: {original_dye_slug})"
                    )

            if dye_slug != original_dye_slug:
                logger.warning(
                    f"Slug collision during Dye migration: {original_dye_slug} -> {dye_slug} "
                    f"(Dye ID: {old_id}, Name: {name})"
                )

            # Old Dye schema doesn't have chemical structure fields
            # Mark all as PROPRIETARY to avoid unique constraint issues
            # (Can be updated later with actual chemical data)

            # Create Dye container (without fluorescence properties)
            dye = Dye.objects.create(
                created=created,
                modified=modified,
                name=name,
                slug=dye_slug,
                inchikey="",  # No chemical data in old schema
                structural_status="PROPRIETARY",  # Safe default for old dyes


                # Note: Other fields like smiles, inchi, etc. can be added later
                # if they exist in the old schema
            )

            # Create Fluorophore for this DyeState
            # Ensure unique fluorophore slug
            fluorophore_slug = f"{dye_slug}-default"
            original_fluor_slug = fluorophore_slug
            base_fluor_slug = fluorophore_slug
            counter = 1
            while Fluorophore.objects.filter(slug=fluorophore_slug).exists():
                fluorophore_slug = f"{base_fluor_slug}-{counter}"
                counter += 1
                if counter > 100:
                    raise ValueError(
                        f"Could not generate unique fluorophore slug for Dye {old_id} "
                        f"after 100 attempts (original: {original_fluor_slug})"
                    )

            if fluorophore_slug != original_fluor_slug:
                logger.warning(
                    f"Fluorophore slug collision during Dye migration: {original_fluor_slug} -> {fluorophore_slug} "
                    f"(Dye ID: {old_id}, Name: {name})"
                )

            # Don't pass emhex/exhex - the save() method will compute them from wavelengths
            fluorophore = Fluorophore.objects.create(
                created=created,
                modified=modified,
                label=name,
                slug=fluorophore_slug,
                entity_type="dye",
                ex_max=ex_max,
                em_max=em_max,
                ext_coeff=ext_coeff,
                qy=qy,
                brightness=brightness,
                lifetime=lifetime,
                pka=pka,
                twop_ex_max=twop_ex_max,
                twop_peakGM=twop_peakgm,
                twop_qy=twop_qy,
                is_dark=is_dark,
                
                
            )

            # Create DyeState (one per old Dye)
            # Use raw SQL to avoid Django MTI trying to update parent with empty values
            cursor.execute("""
                INSERT INTO proteins_dyestate (fluorophore_ptr_id, dye_id, name, solvent, ph, environment, is_reference)
                VALUES (%s, %s, %s, %s, %s, %s, %s)
            """, [fluorophore.pk, dye.pk, "default", "PBS", 7.4, "FREE", True])

            # Create FluorescenceMeasurement from old Dye data
            if any([ex_max, em_max, qy, ext_coeff, lifetime, pka]):
                # For dyes, we don't have a primary_reference concept in old schema
                # We'll leave reference as None for now
                FluorescenceMeasurement.objects.create(
                    fluorophore=fluorophore,
                    reference_id=None,
                    ex_max=ex_max,
                    em_max=em_max,
                    ext_coeff=ext_coeff,
                    qy=qy,
                    brightness=brightness,
                    lifetime=lifetime,
                    pka=pka,
                    twop_ex_max=twop_ex_max,
                    twop_peakGM=twop_peakgm,
                    twop_qy=twop_qy,
                    is_dark=is_dark,
                    is_trusted=True,
                    
                    
                )

    print(f"Migrated {Dye.objects.count()} Dye records to Dye containers")
    print(f"Created {DyeState.objects.count()} DyeState records")


def update_spectrum_ownership(apps, schema_editor):
    """Update Spectrum foreign keys to point to new Fluorophore records."""
    Spectrum = apps.get_model("proteins", "Spectrum")
    Fluorophore = apps.get_model("proteins", "Fluorophore")

    # We need to map old State/Dye IDs to new Fluorophore IDs
    # This is tricky because we need to query the old tables

    with schema_editor.connection.cursor() as cursor:
        # Update spectra that were owned by States
        cursor.execute("""
            UPDATE proteins_spectrum s
            SET owner_fluor_id = (
                SELECT f.id
                FROM proteins_fluorophore f
                JOIN proteins_state ns ON ns.fluorophore_ptr_id = f.id
                JOIN proteins_state_old os ON os.slug = f.slug
                WHERE os.id = s.owner_state_id
            )
            WHERE s.owner_state_id IS NOT NULL
        """)

        state_count = cursor.rowcount
        print(f"Updated {state_count} spectra owned by States")

        # Update spectra that were owned by Dyes
        # Note: DyeState slug is "{dye_slug}-default", so we need to match carefully
        cursor.execute("""
            UPDATE proteins_spectrum s
            SET owner_fluor_id = (
                SELECT f.id
                FROM proteins_fluorophore f
                JOIN proteins_dyestate ds ON ds.fluorophore_ptr_id = f.id
                JOIN proteins_dye d ON d.id = ds.dye_id
                JOIN proteins_dye_old od ON od.slug = d.slug
                WHERE od.id = s.owner_dye_id
            )
            WHERE s.owner_dye_id IS NOT NULL
        """)

        dye_count = cursor.rowcount
        print(f"Updated {dye_count} spectra owned by Dyes")


def update_ocfluoreff(apps, schema_editor):
    """Update OcFluorEff to use direct FK to Fluorophore."""
    with schema_editor.connection.cursor() as cursor:
        # Update OcFluorEff records that pointed to States via GenericFK
        cursor.execute("""
            UPDATE proteins_ocfluoreff o
            SET fluor_id = (
                SELECT f.id
                FROM proteins_fluorophore f
                JOIN proteins_state s ON s.fluorophore_ptr_id = f.id
                WHERE s.fluorophore_ptr_id IN (
                    SELECT fluorophore_ptr_id
                    FROM proteins_state ps
                    JOIN proteins_state_old os ON os.slug = (
                        SELECT slug FROM proteins_fluorophore WHERE id = ps.fluorophore_ptr_id
                    )
                    WHERE os.id = o.object_id::integer
                )
            )
            WHERE o.content_type_id = (
                SELECT id FROM django_content_type
                WHERE app_label = 'proteins' AND model = 'state'
            )
        """)

        state_count = cursor.rowcount
        print(f"Updated {state_count} OcFluorEff records that pointed to States")

        # Update OcFluorEff records that pointed to Dyes via GenericFK
        cursor.execute("""
            UPDATE proteins_ocfluoreff o
            SET fluor_id = (
                SELECT f.id
                FROM proteins_fluorophore f
                JOIN proteins_dyestate ds ON ds.fluorophore_ptr_id = f.id
                JOIN proteins_dye d ON d.id = ds.dye_id
                JOIN proteins_dye_old od ON od.slug = d.slug
                WHERE od.id = o.object_id::integer
            )
            WHERE o.content_type_id = (
                SELECT id FROM django_content_type
                WHERE app_label = 'proteins' AND model = 'dye'
            )
        """)

        dye_count = cursor.rowcount
        print(f"Updated {dye_count} OcFluorEff records that pointed to Dyes")


def migrate_forward(apps, schema_editor):
    """Run all migration functions."""
    print("Starting data migration from old schema...")
    migrate_state_data(apps, schema_editor)
    migrate_dye_data(apps, schema_editor)
    update_spectrum_ownership(apps, schema_editor)
    update_ocfluoreff(apps, schema_editor)
    print("Data migration complete!")


def migrate_reverse(apps, schema_editor):
    """Reverse migration - not supported.

    This migration performs a one-way data transformation from the old schema
    (separate State and Dye models) to the new schema (MTI-based Fluorophore hierarchy).

    Reversing this migration would require:
    1. Decomposing Fluorophore + State back into old State structure
    2. Decomposing Fluorophore + DyeState back into old Dye structure
    3. Merging FluorescenceMeasurement data back into parent entities
    4. Restoring old spectrum ownership relationships

    This is not safely automatable. If you need to rollback this migration:
    1. Restore from a database backup taken before running this migration
    2. Do NOT attempt to use Django's reverse migration functionality
    3. Estimated restore time: 30-60 minutes depending on database size

    See deployment documentation for rollback procedures.
    """
    raise RuntimeError(
        "This migration cannot be reversed. "
        "Restore from database backup if rollback is needed. "
        "See migration docstring for details."
    )


class Migration(migrations.Migration):
    dependencies = [
        ("proteins", "0059_add_fluorophore_and_new_models"),
    ]

    operations = [
        migrations.RunPython(migrate_forward, migrate_reverse),
    ]
