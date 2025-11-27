# Spectrum Model Evaluation

Executive Summary

The Spectrum model is well-designed for its purpose but has several areas where improvements could yield significant benefits. The current design prioritizes correctness and
flexibility over storage efficiency, which is reasonable for a ~7,500 row table. However, there are some data integrity edge cases and optimization opportunities worth
considering.

---

1. Field Design & Data Integrity

Current Design ✅ Strengths

| Field                 | Assessment                                                                        |
|-----------------------|-----------------------------------------------------------------------------------|
| category / subtype    | Good choice with db_index. Enum-like choices with validation                      |
| status (StatusModel)  | Clean approval workflow with proper indexing                                      |
| Polymorphic ownership | Explicit FKs avoid GenericForeignKey complexity (good decision citing Luke Plant) |
| reference FK          | Proper nullable FK for attribution                                                |

⚠️ Integrity Concerns Found

1. Two-photon normalization has edge cases:
ID 9902: range [0.0001, 1.1310]  ← max > 1.0
ID 192:  range [0.0006, 10.4490] ← severely out of range
ID 193:  range [0.0000, 4.4754]
ID 160:  range [0.0167, 3.5225]
1. The norm2P() function uses peak detection with order=100, which can fail for noisy or unusual spectra. There's no post-validation that values are actually in [0,1].
2. Missing database-level constraints: The single-owner validation happens in Python (save()) but not at the DB level:
if sum(bool(x) for x in self.owner_set) > 1:
    raise ValidationError("Spectrum must have only one owner!")
2. A CHECK constraint would prevent data corruption from raw SQL or bulk operations.
3. Precision stored exceeds need: Max decimal places = 20 (float64 artifacts), but scientific precision needs only 4-5 decimals. Storing 4 decimal places is sufficient for
0.01% precision.

Recommendations

-- Add database-level constraints
ALTER TABLE proteins_spectrum ADD CONSTRAINT spectrum_single_owner CHECK (
    (CASE WHEN owner_fluor_id IS NOT NULL THEN 1 ELSE 0 END +
    CASE WHEN owner_filter_id IS NOT NULL THEN 1 ELSE 0 END +
    CASE WHEN owner_light_id IS NOT NULL THEN 1 ELSE 0 END +
    CASE WHEN owner_camera_id IS NOT NULL THEN 1 ELSE 0 END) = 1
);

-- Add min_wave/max_wave as computed stored columns (PostgreSQL 12+)
ALTER TABLE proteins_spectrum
ADD COLUMN min_wave SMALLINT GENERATED ALWAYS AS ([data[1]](1)::smallint) STORED,
ADD COLUMN max_wave SMALLINT GENERATED ALWAYS AS ((data[array_length(data,1)])[1]::smallint) STORED;

---

2. Read-Heavy Query Optimization

Current Performance

| Query Type                | Time                |
|---------------------------|---------------------|
| Metadata-only (1000 rows) | 3.3ms               |
| Full data (1000 rows)     | 339ms ← 100× slower |

The data field is the bottleneck. PostgreSQL must deserialize the entire 2D float array even if you only need metadata.

Current Optimizations ✅

- get_spectra_list() uses .only() and .values() to avoid loading data
- GraphQL uses SpectrumInfo type for metadata-only queries
- Individual spectra cached with 24-hour TTL
- Bulk spectra list cached indefinitely with signal-based invalidation

Missing Optimizations ⚠️

1. No TOAST threshold tuning: PostgreSQL stores arrays inline up to ~2KB, then TOASTs. With avg 11KB per spectrum, all are TOASTed, requiring separate I/O.
2. No covering indexes for common queries: The spectra viewer fetches metadata without data, but indexes don't cover all needed columns.
3. Computed columns not leveraged: min_wave, max_wave, peak_wave are computed on every access via Python properties.

Recommendations

-- Covering index for metadata queries (avoids table access)
CREATE INDEX spectrum_metadata_covering_idx ON proteins_spectrum (
    status, category, subtype, owner_fluor_id, owner_filter_id, owner_light_id, owner_camera_id
);

-- Store computed values (min_wave is data[1][1], max_wave is data[-1][1])
ALTER TABLE proteins_spectrum
ADD COLUMN min_wave SMALLINT,
ADD COLUMN max_wave SMALLINT,
ADD COLUMN peak_wave SMALLINT;

---

3. PostgreSQL-Specific Improvements

Current Usage

- ArrayField (PostgreSQL-specific) ✅
- TrigramSimilarity for fuzzy search ✅
- Proper indexing strategy ✅

Underutilized PostgreSQL Features

1. BYTEA for compact storage: Since all wavelengths are integers at 1nm steps, you only need to store Y-values:
-- Current: float8[][] at ~17 bytes per point
-- Alternative: bytea with uint16 Y-values at 2 bytes per point
-- 88% storage reduction
2. Generated columns (PostgreSQL 12+):
min_wave SMALLINT GENERATED ALWAYS AS ([data[1]](1)::smallint) STORED
3. Partial indexes for common filters:
CREATE INDEX spectrum_approved_fluor_idx ON proteins_spectrum (owner_fluor_id)
WHERE status = 'approved' AND owner_fluor_id IS NOT NULL;
4. BRIN index on created/modified for time-range queries (if needed):
CREATE INDEX spectrum_created_brin ON proteins_spectrum USING BRIN (created);

---

4. Data Normalization Analysis

Current Normalization Pipeline

Input Data → Validate → Interpolate to 1nm → Normalize to peak=1 → Store

Should You Be Normalizing?

Arguments FOR normalization (current approach):

- Consistent 1nm grid enables direct comparisons/overlays without runtime interpolation
- Peak normalization to 1.0 is scientifically meaningful (relative intensities)
- Reduces client-side complexity for visualization
- Enables efficient database operations (fixed array sizes per wavelength range)

Arguments AGAINST:

- Loss of original data (can't recover absolute intensities without EC/QY)
- Some two-photon spectra have peaks >1.0 (normalization bugs)
- Interpolation can introduce artifacts, especially for sharp features
- 1nm resolution may be overkill for broad spectra (could use 2nm for storage savings)

Verdict: Keep normalizing, but improve robustness

The 1nm normalization is scientifically sound and enables the core use case (spectral overlap visualization). However:

1. Store original peak value (you partially do this for 2P with _peakval2p, but it's not persisted)
2. Add post-normalization validation to catch the >1.0 edge cases
3. Consider storing original data separately for provenance (in a reference field or separate table)

---

5. Storage Optimization Options

Current Storage Analysis

| Metric                  | Value                          |
|-------------------------|--------------------------------|
| Total spectra           | 7,478                          |
| Table size              | 3 MB (data) + 1.5 MB (indexes) |
| Avg JSON representation | 11,265 bytes                   |
| Avg data points         | 684 points                     |

Alternative Storage Approaches

| Approach                      | Size/Spectrum | Savings | Complexity |
|-------------------------------|---------------|---------|------------|
| Current (float8[][])          | ~11KB         | —       | Low        |
| Y-only + min_wave             | ~6KB          | 45%     | Low        |
| Binary float64                | ~4.5KB        | 60%     | Medium     |
| Binary uint16 (0-10000 → 0-1) | ~1.1KB        | 90%     | Medium     |
| Delta-encoded + zstd          | ~0.5KB        | 95%     | High       |

Recommended Approach: Binary Y-values with metadata

Since wavelengths are always integers at 1nm steps:

# New schema concept

class Spectrum(models.Model):
    min_wave = models.SmallIntegerField()  # e.g., 350
    # Y-values stored as bytea: 2 bytes per point, values 0-10000 → 0.0-1.0
    y_data = models.BinaryField()

    @property
    def data(self):
        """Reconstruct [[wave, val], ...] format for backward compatibility."""
        import struct
        y_vals = struct.unpack(f'{len(self.y_data)//2}H', self.y_data)
        return [[self.min_wave + i, v / 10000] for i, v in enumerate(y_vals)]

Benefits:

- 90% storage reduction (3MB → 300KB for spectra table)
- Faster I/O (less data to read/deserialize)
- Maintains backward compatibility via property
- 4 decimal precision (0.0001) is sufficient for normalized spectra

Trade-offs:

- Migration complexity
- Slightly more CPU for encoding/decoding (negligible)
- Raw SQL queries become harder (but you use ORM)

---

6. Domain-Specific Recommendations

As someone who works with spectral data:

Keep Doing ✅

1. 1nm interpolation - Standard for spectral databases, enables overlap calculations
2. Peak normalization - Essential for comparing relative spectral shapes
3. Separate EX/ABS subtypes - Scientifically distinct (EX is functionally relevant, ABS is raw measurement)
4. pH/solvent metadata - Critical for rigorous spectral comparisons

Consider Adding

1. Instrument/conditions metadata - Bandwidth, integration time, detector type
2. Uncertainty/noise estimate - Even a simple SNR value helps quality assessment
3. Original data hash - For verifying data hasn't been corrupted
4. Measurement temperature - Spectra are temperature-dependent

Fix These Issues

1. Two-photon normalization bugs - 4 spectra have values >1.0
2. Absorption vs Excitation distinction - Some proteins have both, ensure they're not confused
3. Filter transmission validation - Ensure filter spectra are in 0-1 range (not 0-100%)

---
Summary of Recommendations

Quick Wins (Low Effort, High Value)

1. Add database CHECK constraint for single-owner validation
2. Add post-normalization validation in clean() to ensure values ≤ 1.0
3. Store min_wave/max_wave as computed columns (avoid property calls)

Medium-Term (Moderate Effort)

1. Create covering index for metadata-only queries
2. Add partial indexes for status='approved' queries
3. Fix the 4 broken two-photon spectra

Long-Term (If Performance Becomes an Issue)

1. Migrate to binary Y-value storage (90% size reduction)
2. Consider separate table for raw/original data
3. Implement BRIN indexes if time-based queries become common
