from __future__ import annotations

from django.core.files.base import ContentFile

from fpbase.storage import ViteManifestStaticFilesStorage


def test_vite_assets_keep_their_names(tmp_path):
    storage = ViteManifestStaticFilesStorage(location=tmp_path)
    content = ContentFile(b"x")

    assert storage.hashed_name("assets/icons-D265L9a5.js", content) == "assets/icons-D265L9a5.js"
    # everything else still gets Django's content hash
    assert storage.hashed_name("images/logo.png", content) != "images/logo.png"
