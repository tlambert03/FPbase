# Typing Migration Notes

This document tracks questions, uncertainties, and decisions made during the migration to stricter type checking with `ty` and `django-types`.

## Summary

Successfully migrated FPbase to use `ty` type checker with `django-types`. Started with 57 errors (after excluding migrations and tests), fixed actual bugs, and added appropriate type ignores for third-party library stubs issues.

## Decisions Made

### Configuration
- Excluded migrations and tests from type checking (configured in `pyproject.toml`)
- Set line length to 119 chars (project standard)
- Ignored `possibly-unbound-attribute` rule globally (Django patterns trigger this)

### Actual Code Fixes Made
1. **fpseq/mutations.py**: Added None check for `muts` parameter to prevent iteration over None
2. **fpseq/skbio_protein.py**: Replaced deprecated `np.in1d` with `np.isin` (3 occurrences)
3. **references/models.py**: Added None check for `authors` before iteration
4. **proteins/views/protein.py**: Changed `gbseqs.get()` to direct dict access `gbseqs[]` since key existence was already checked
5. **proteins/util/helpers.py**: Added None check for `slug_dict` cache value
6. **proteins/util/helpers.py**: Converted matplotlib position list to tuple with ignore comment
7. **proteins/util/spectra.py**: Converted `range` to `np.array` for scipy interpolation

### Type Ignores Added

#### Django Stubs Issues
- **Form Meta inheritance** (2): `UserChangeForm.Meta` and `UserCreationForm.Meta` not in stubs
- **Django forms.widgets** (10+): Module path not recognized in django-types
- **DRF Throttled.wait** (1): Attribute exists at runtime but not in stubs
- **BaseModelFormSet.form** (2): Dynamic formset attribute not typed
- **models.ObjectDoesNotExist** (1): Exception class path issue
- **model.STATUS** (2): Dynamic enum attribute on models
- **Form widget.render()** (1): Signature mismatch in stubs
- **ProteinFilter.recs** (1): Dynamic filter attribute

#### Third-Party Library Issues
- **django_filters.rest_framework** (2): Possibly-unbound import warnings
- **graphene_django.filter** (1): Possibly-unbound import
- **graphene.relay** (1): Module attribute not recognized
- **graphene.Enum** (1): Call signature issue
- **BioPython** (2): Missing stubs for `Data.CodonTable` and `Seq` constructor

#### NumPy/SciPy Issues
- **NumPy array iteration** (2): Type checker can't verify array is iterable
- **callable type hint** (1): Built-in `callable` not allowed in type expressions

## Errors Encountered

### Initial State
- **57 errors** after excluding migrations and tests
- Breakdown: ~30% actual bugs, ~70% third-party stubs issues

### Final State
- **0 errors** - All checks passing
- **All 174 tests passing** (1 xfailed, unrelated to typing changes)

## Notes for Future

### Potential Improvements
1. Consider contributing fixes to django-types for missing attributes
2. Investigate if Factory Boy has better type support we could use
3. Some ignores could be avoided with better type annotations on our models
4. Consider using `from __future__ import annotations` more widely

### Things to Watch
- When upgrading django-types, some ignores may no longer be needed
- NumPy type stubs continue to improve, array iteration issues may resolve
- ty is in alpha - expect rapid improvements and potential new features
