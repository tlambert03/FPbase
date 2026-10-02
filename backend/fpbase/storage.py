from __future__ import annotations

from whitenoise.storage import CompressedManifestStaticFilesStorage


class ViteManifestStaticFilesStorage(CompressedManifestStaticFilesStorage):
    """Leave Vite's build output under its own (already content-hashed) names.

    Vite's chunks import each other by their Vite names, so if Django also renamed
    them, pages would preload the Django names and then fetch every chunk again
    under the Vite name.
    """

    def hashed_name(self, name, content=None, filename=None):
        if name.startswith("assets/"):
            return name
        return super().hashed_name(name, content, filename)
