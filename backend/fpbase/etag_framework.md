# ETag Framework Architecture

## Goal
Provide reusable HTTP ETag support for GraphQL, REST API, and HTML views to enable conditional requests (304 Not Modified responses).

## Design Principles
1. **Simple**: Easy to apply to any view with minimal code
2. **Automatic**: Signal-based cache invalidation, no manual calls
3. **Flexible**: Works for GraphQL, DRF, and Django template views
4. **Efficient**: Uses Redis cache for version tracking

## Components

### 1. Model Version Tracking (`fpbase/cache_utils.py`)

**Purpose**: Track when models change using cache-based versioning.

**API**:
```python
get_model_version(*model_classes) -> str
    """Get combined version hash for multiple models."""

invalidate_model_version(model_class)
    """Bump version when model changes (called by signals)."""
```

**Storage**: Redis cache keys: `model_version:app_label.ModelName`

**Signals**: Auto-registered for Spectrum, Protein, State, OpticalConfig, etc.

### 2. ETag Generation (`fpbase/etag_utils.py`)

**Purpose**: Generate ETags from model versions or response content.

**API**:
```python
generate_version_etag(*model_classes) -> str
    """Generate weak ETag from model versions."""
    # Example: W/"abc123def456"

generate_content_etag(content: bytes | str) -> str
    """Generate strong ETag from response content."""
    # Example: "d41d8cd98f00b204e9800998ecf8427e"

parse_etag_header(header_value: str) -> list[str]
    """Parse If-None-Match header into list of ETags."""
```

**ETag Types**:
- **Weak** (`W/"..."`) for version-based (data hasn't changed)
- **Strong** (`"..."`) for content-based (byte-for-byte identical)

### 3. Conditional Request Handling

#### For DRF Views (`fpbase/api_mixins.py`)
```python
class ETagMixin:
    """Add ETag support to DRF APIView classes."""
    etag_models = []  # Define in view: [Protein, State]

    def finalize_response(self, request, response, *args, **kwargs):
        # Add ETag header
        # Check If-None-Match
        # Return 304 if matched
```

#### For GraphQL (`fpbase/views.py`)
```python
class RateLimitedGraphQLView(GraphQLView):
    def dispatch(self, request, *args, **kwargs):
        # Check If-None-Match for Spectrum + OpticalConfig
        # Add ETag to response
```

## Usage Examples

### REST API
```python
from fpbase.api_mixins import ETagMixin
from proteins.models import Protein, State

class ProteinListAPIView(ETagMixin, ListAPIView):
    etag_models = [Protein, State]  # Auto-handled!
    # ... existing code ...
```

### GraphQL
```python
# Already integrated in RateLimitedGraphQLView
# Tracks Spectrum + OpticalConfig versions
```

## Cache Invalidation Flow

```
1. User saves Spectrum instance
2. post_save signal fires
3. invalidate_model_version(Spectrum) called
4. Cache key "model_version:proteins.Spectrum" updated to new timestamp
5. Next request generates new ETag (old client ETags won't match)
6. Client receives new data with new ETag
```

## Testing Strategy

1. **Unit tests**: Version tracking, ETag generation
2. **Integration tests**: Signal handlers, cache invalidation
3. **View tests**: 200 vs 304 responses, header handling
4. **E2E tests**: Real browser caching behavior

## Performance Impact

- **Cache reads**: 1-2 per request (model versions)
- **Cache writes**: Only on model changes (rare)
- **Response size savings**: 95%+ for cached responses (304 = ~300 bytes)
- **Server CPU savings**: No serialization for 304 responses
