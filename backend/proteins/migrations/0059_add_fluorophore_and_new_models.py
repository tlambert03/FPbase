# Generated manually for schema overhaul

from django.core.validators import MaxValueValidator, MinValueValidator
from django.contrib.postgres.fields import ArrayField
from django.db import migrations, models
import django.db.models.deletion
import model_utils.fields
import django.utils.timezone


class Migration(migrations.Migration):
    dependencies = [
        ("proteins", "0058_snapgeneplasmid_protein_snapgene_plasmids"),
        ("references", "0001_initial"),
    ]

    operations = [
        migrations.SeparateDatabaseAndState(
            database_operations=[
                # In the database: copy old tables to _old versions, then drop originals
                # This avoids index naming conflicts
                migrations.RunSQL(
                    sql="""
                        CREATE TABLE proteins_state_old AS SELECT * FROM proteins_state;
                        DROP TABLE proteins_state CASCADE;
                    """,
                    reverse_sql="DROP TABLE IF EXISTS proteins_state_old CASCADE;",
                ),
                migrations.RunSQL(
                    sql="""
                        CREATE TABLE proteins_dye_old AS SELECT * FROM proteins_dye;
                        DROP TABLE proteins_dye CASCADE;
                    """,
                    reverse_sql="DROP TABLE IF EXISTS proteins_dye_old CASCADE;",
                ),
            ],
            state_operations=[
                # In Django's state: remove old models
                migrations.DeleteModel(name="State"),
                migrations.DeleteModel(name="Dye"),
            ],
        ),

        # Now create all new models
        migrations.CreateModel(
            name="Fluorophore",
            fields=[
                ("id", models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name="ID")),
                ("created", model_utils.fields.AutoCreatedField(default=django.utils.timezone.now, editable=False, verbose_name="created")),
                ("modified", model_utils.fields.AutoLastModifiedField(default=django.utils.timezone.now, editable=False, verbose_name="modified")),
                # Fluorescence data fields
                ("ex_max", models.PositiveSmallIntegerField(blank=True, null=True, validators=[MinValueValidator(300), MaxValueValidator(900)], db_index=True, help_text="Excitation maximum (nm)")),
                ("em_max", models.PositiveSmallIntegerField(blank=True, null=True, validators=[MinValueValidator(300), MaxValueValidator(1000)], db_index=True, help_text="Emission maximum (nm)")),
                ("emhex", models.CharField(max_length=7, blank=True)),
                ("exhex", models.CharField(max_length=7, blank=True)),
                ("ext_coeff", models.IntegerField(blank=True, null=True, verbose_name="Extinction Coefficient (M-1 cm-1)", validators=[MinValueValidator(0), MaxValueValidator(300000)])),
                ("qy", models.FloatField(null=True, blank=True, verbose_name="Quantum Yield", validators=[MinValueValidator(0), MaxValueValidator(1)])),
                ("brightness", models.FloatField(null=True, blank=True, editable=False)),
                ("lifetime", models.FloatField(null=True, blank=True, help_text="Lifetime (ns)", validators=[MinValueValidator(0), MaxValueValidator(20)])),
                ("pka", models.FloatField(null=True, blank=True, verbose_name="pKa", validators=[MinValueValidator(2), MaxValueValidator(12)])),
                ("twop_ex_max", models.PositiveSmallIntegerField(blank=True, null=True, verbose_name="Peak 2P excitation", validators=[MinValueValidator(700), MaxValueValidator(1600)], db_index=True)),
                ("twop_peakGM", models.FloatField(null=True, blank=True, verbose_name="Peak 2P cross-section of S0->S1 (GM)", validators=[MinValueValidator(0), MaxValueValidator(200)])),
                ("twop_qy", models.FloatField(null=True, blank=True, verbose_name="2P Quantum Yield", validators=[MinValueValidator(0), MaxValueValidator(1)])),
                ("is_dark", models.BooleanField(default=False, verbose_name="Dark State", help_text="This state does not fluorescence")),
                # Identity fields
                ("label", models.CharField(max_length=255, db_index=True)),
                ("slug", models.SlugField(max_length=100, unique=True)),  # Increased from default 50 to accommodate long dye names
                ("entity_type", models.CharField(max_length=10, choices=[("protein", "Protein"), ("dye", "Dye")], db_index=True)),
                # Lineage tracking
                ("source_map", models.JSONField(default=dict, blank=True)),
                # Author tracking (from Authorable mixin)
                ("created_by_id", models.IntegerField(null=True, blank=True)),
                ("updated_by_id", models.IntegerField(null=True, blank=True)),
            ],
            options={
                "indexes": [
                    models.Index(fields=["ex_max"], name="fluorophore_ex_max_idx"),
                    models.Index(fields=["em_max"], name="fluorophore_em_max_idx"),
                    models.Index(fields=["label", "entity_type"], name="fluorophore_label_type_idx"),
                    models.Index(fields=["entity_type", "is_dark"], name="fluorophore_type_dark_idx"),
                ],
            },
        ),

        migrations.CreateModel(
            name="FluorescenceMeasurement",
            fields=[
                ("id", models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name="ID")),
                ("created", model_utils.fields.AutoCreatedField(default=django.utils.timezone.now, editable=False, verbose_name="created")),
                ("modified", model_utils.fields.AutoLastModifiedField(default=django.utils.timezone.now, editable=False, verbose_name="modified")),
                # Fluorescence data fields (same as Fluorophore)
                ("ex_max", models.PositiveSmallIntegerField(blank=True, null=True, validators=[MinValueValidator(300), MaxValueValidator(900)], db_index=True, help_text="Excitation maximum (nm)")),
                ("em_max", models.PositiveSmallIntegerField(blank=True, null=True, validators=[MinValueValidator(300), MaxValueValidator(1000)], db_index=True, help_text="Emission maximum (nm)")),
                ("emhex", models.CharField(max_length=7, blank=True)),
                ("exhex", models.CharField(max_length=7, blank=True)),
                ("ext_coeff", models.IntegerField(blank=True, null=True, verbose_name="Extinction Coefficient (M-1 cm-1)", validators=[MinValueValidator(0), MaxValueValidator(300000)])),
                ("qy", models.FloatField(null=True, blank=True, verbose_name="Quantum Yield", validators=[MinValueValidator(0), MaxValueValidator(1)])),
                ("brightness", models.FloatField(null=True, blank=True, editable=False)),
                ("lifetime", models.FloatField(null=True, blank=True, help_text="Lifetime (ns)", validators=[MinValueValidator(0), MaxValueValidator(20)])),
                ("pka", models.FloatField(null=True, blank=True, verbose_name="pKa", validators=[MinValueValidator(2), MaxValueValidator(12)])),
                ("twop_ex_max", models.PositiveSmallIntegerField(blank=True, null=True, verbose_name="Peak 2P excitation", validators=[MinValueValidator(700), MaxValueValidator(1600)], db_index=True)),
                ("twop_peakGM", models.FloatField(null=True, blank=True, verbose_name="Peak 2P cross-section of S0->S1 (GM)", validators=[MinValueValidator(0), MaxValueValidator(200)])),
                ("twop_qy", models.FloatField(null=True, blank=True, verbose_name="2P Quantum Yield", validators=[MinValueValidator(0), MaxValueValidator(1)])),
                ("is_dark", models.BooleanField(default=False, verbose_name="Dark State", help_text="This state does not fluorescence")),
                # Measurement-specific fields
                ("date_measured", models.DateField(null=True, blank=True)),
                ("conditions", models.TextField(blank=True, help_text="pH, solvent, temp, etc.")),
                ("is_trusted", models.BooleanField(default=False, help_text="If True, this measurement overrides others.")),
                # Foreign keys
                ("fluorophore", models.ForeignKey("Fluorophore", related_name="measurements", on_delete=django.db.models.deletion.CASCADE)),
                ("reference", models.ForeignKey("references.Reference", on_delete=django.db.models.deletion.CASCADE, null=True, blank=True)),
                # Author tracking fields
                ("created_by_id", models.IntegerField(null=True, blank=True)),
                ("updated_by_id", models.IntegerField(null=True, blank=True)),
            ],
            options={
                "abstract": False,
            },
        ),

        migrations.CreateModel(
            name="Dye",
            fields=[
                ("id", models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name="ID")),
                ("created", model_utils.fields.AutoCreatedField(default=django.utils.timezone.now, editable=False, verbose_name="created")),
                ("modified", model_utils.fields.AutoLastModifiedField(default=django.utils.timezone.now, editable=False, verbose_name="modified")),
                ("name", models.CharField(max_length=255, db_index=True)),
                ("slug", models.SlugField(max_length=100, unique=True)),  # Increased from default 50 to accommodate long names
                ("synonyms", ArrayField(models.CharField(max_length=255), blank=True, default=list)),
                ("structural_status", models.CharField(max_length=20, choices=[("DEFINED", "Defined Structure"), ("PROPRIETARY", "Proprietary / Unknown Structure")], default="DEFINED")),
                ("canonical_smiles", models.TextField(blank=True)),
                ("inchi", models.TextField(blank=True)),
                ("inchikey", models.CharField(max_length=27, blank=True, db_index=True)),
                ("molblock", models.TextField(blank=True, help_text="V3000 Molfile for precise rendering")),
                ("parent_mixture", models.ForeignKey("self", on_delete=django.db.models.deletion.SET_NULL, null=True, blank=True, related_name="isomers")),
                ("chemical_class", models.CharField(max_length=100, blank=True, db_index=True)),
                ("equilibrium_constant_klz", models.FloatField(null=True, blank=True, help_text="Equilibrium constant between non-fluorescent lactone and fluorescent zwitterion.")),
                # Product mixin fields
                ("manufacturer", models.CharField(max_length=128, blank=True)),
                ("part", models.CharField(max_length=128, blank=True)),
                ("url", models.URLField(blank=True)),
                # Author tracking fields
                ("created_by_id", models.IntegerField(null=True, blank=True)),
                ("updated_by_id", models.IntegerField(null=True, blank=True)),
            ],
            options={
                "constraints": [
                    models.UniqueConstraint(
                        fields=["inchikey"],
                        name="unique_defined_molecule",
                        condition=models.Q(structural_status="DEFINED"),
                    )
                ],
            },
        ),

        migrations.CreateModel(
            name="DyeState",
            fields=[
                ("fluorophore_ptr", models.OneToOneField(auto_created=True, on_delete=django.db.models.deletion.CASCADE, parent_link=True, primary_key=True, serialize=False, to="proteins.fluorophore")),
                ("name", models.CharField(max_length=255, help_text="e.g., 'Bound to DNA' or 'In Methanol'")),
                ("solvent", models.CharField(max_length=100, default="PBS")),
                ("ph", models.FloatField(default=7.4)),
                ("environment", models.CharField(max_length=20, choices=[], default="FREE")),
                ("is_reference", models.BooleanField(default=False, help_text="If True, this is the default state shown on the dye summary card.")),
                ("dye", models.ForeignKey("Dye", on_delete=django.db.models.deletion.CASCADE, related_name="states")),
            ],
            options={
                "abstract": False,
            },
            bases=("proteins.fluorophore",),
        ),

        # Re-create State model as MTI child of Fluorophore
        migrations.CreateModel(
            name="State",
            fields=[
                ("fluorophore_ptr", models.OneToOneField(auto_created=True, on_delete=django.db.models.deletion.CASCADE, parent_link=True, primary_key=True, serialize=False, to="proteins.fluorophore")),
                ("name", models.CharField(max_length=64, default="default")),
                ("protein", models.ForeignKey("Protein", related_name="states", help_text="The protein to which this state belongs", on_delete=django.db.models.deletion.CASCADE)),
                ("maturation", models.FloatField(null=True, blank=True, help_text="Maturation time (min)", validators=[MinValueValidator(0), MaxValueValidator(1600)])),
            ],
            options={
                "abstract": False,
            },
            bases=("proteins.fluorophore",),
        ),

        # Add owner_fluor to Spectrum (nullable for now)
        migrations.AddField(
            model_name="spectrum",
            name="owner_fluor",
            field=models.ForeignKey(
                blank=True,
                null=True,
                on_delete=django.db.models.deletion.CASCADE,
                related_name="spectra",
                to="proteins.fluorophore",
            ),
        ),

        # Add fluor FK to OcFluorEff (nullable for now)
        migrations.AddField(
            model_name="ocfluoreff",
            name="fluor",
            field=models.ForeignKey(
                blank=True,
                null=True,
                on_delete=django.db.models.deletion.CASCADE,
                related_name="oc_effs",
                to="proteins.fluorophore",
            ),
        ),
    ]
