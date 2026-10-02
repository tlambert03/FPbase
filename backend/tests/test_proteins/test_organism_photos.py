from __future__ import annotations

import json
from pathlib import Path

import pytest
from django.core.exceptions import ValidationError

from proteins.factories import OrganismFactory
from proteins.models import OrganismPhoto
from proteins.models.organism import COMMONS_THUMB_WIDTHS

MIGRATION_DATA = (
    Path(__file__).parents[2] / "proteins" / "migrations" / "0066_organism_photos.json"
)


def _photo(organism, **kwargs) -> OrganismPhoto:
    defaults = {
        "commons_file": "Acropora digitifera, Erub.JPG",
        "width": 3456,
        "height": 2592,
        "author": "Kerryn Johns",
        "license": "CC BY 3.0",
        "license_url": "https://creativecommons.org/licenses/by/3.0/",
    }
    return OrganismPhoto.objects.create(organism=organism, **{**defaults, **kwargs})


@pytest.mark.django_db
def test_commons_urls():
    photo = _photo(OrganismFactory(id=70779, scientific_name="Acropora digitifera"))
    # the hashed path Commons itself serves for this file
    base = "https://upload.wikimedia.org/wikipedia/commons"
    name = "Acropora_digitifera%2C_Erub.JPG"
    assert photo.url() == f"{base}/d/db/{name}"
    assert photo.url(500) == f"{base}/thumb/d/db/{name}/500px-{name}"
    assert photo.credit_source_url == (
        "https://commons.wikimedia.org/wiki/File:Acropora_digitifera,_Erub.JPG"
    )
    assert photo.srcset().count("px-") == len(COMMONS_THUMB_WIDTHS)
    assert photo.size_at(960) == (960, 720)


@pytest.mark.django_db
def test_small_commons_image_is_never_upscaled():
    photo = _photo(OrganismFactory(id=3702), width=200, height=600)
    assert photo.src == photo.url()  # the original, not a 960px thumbnail
    assert photo.size_at(960) == (200, 600)
    assert photo.srcset() == f"{photo.url()} 200w"


@pytest.mark.django_db
def test_photo_needs_a_source():
    photo = OrganismPhoto(organism=OrganismFactory(id=6100), author="x", license="CC0")
    with pytest.raises(ValidationError):
        photo.clean()


@pytest.mark.django_db
def test_organism_detail_shows_photo_and_credit(client):
    organism = OrganismFactory(id=86600, scientific_name="Discosoma sp.")
    photo = _photo(organism, pictured="Discosoma nummiforme", pictured_note="same genus")

    html = client.get(organism.get_absolute_url()).content.decode()
    assert photo.src in html
    assert f'<a href="{photo.credit_source_url}"' in html
    assert ">Acropora digitifera, Erub</a>” by" in html  # the work's title
    assert "Kerryn Johns" in html
    assert photo.license_url in html
    assert photo.credit_source_url in html
    assert "Pictured: <em>Discosoma nummiforme</em> (same genus)" in html


@pytest.mark.django_db
def test_organism_pages_without_photo(client):
    organism = OrganismFactory(id=1076, scientific_name="Rhodopseudomonas palustris")
    response = client.get(organism.get_absolute_url())
    assert response.status_code == 200
    assert "organism-photo" not in response.content.decode()

    response = client.get("/organisms/")
    assert response.status_code == 200
    assert "organism-card-placeholder" in response.content.decode()


@pytest.mark.django_db
def test_organism_list_thumbnails(client):
    """The grid links each thumbnail to its organism page, where the full credit is."""
    organism = OrganismFactory(id=70779, scientific_name="Acropora digitifera")
    photo = _photo(organism)
    html = client.get("/organisms/").content.decode()
    assert photo.thumb_src in html
    assert 'title="Photo: Kerryn Johns · CC BY 3.0"' in html
    assert organism.get_absolute_url() in html


def test_migration_data_is_complete():
    photos = json.loads(MIGRATION_DATA.read_text())
    ids = [p["organism_id"] for p in photos]
    assert len(ids) == len(set(ids))
    for p in photos:
        assert p["commons_file"] and p["width"] and p["height"] and p["author"]
        assert p["license"]
        # everything but public domain / CC0 needs a link to its license
        if p["license"] not in ("Public domain", "CC0"):
            assert p["license_url"].startswith("https://creativecommons.org/")
