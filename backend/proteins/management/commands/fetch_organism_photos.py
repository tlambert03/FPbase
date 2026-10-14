from __future__ import annotations

import time
from io import BytesIO
from pathlib import Path

import httpx
from django.core.files.base import ContentFile
from django.core.management.base import BaseCommand
from PIL import Image

from proteins.models import OrganismPhoto

# Wikimedia's robot policy asks for a descriptive User-Agent with a way to reach us
USER_AGENT = "FPbaseBot/1.0 (https://www.fpbase.org/contact/)"
IMAGE_WIDTH = 960  # one of Wikimedia's standard thumbnail sizes
THUMB_WIDTH = 330


class Command(BaseCommand):
    help = "Copy Wikimedia Commons organism photos into our own storage, with thumbnails."

    def handle(self, **options):
        photos = OrganismPhoto.objects.select_related("organism")
        with httpx.Client(headers={"User-Agent": USER_AGENT}, http2=True, timeout=60) as client:
            for photo in photos:
                if photo.image and photo.thumbnail:
                    continue
                if not photo.image:
                    if not photo.commons_file:
                        continue
                    self.stdout.write(f"Downloading {photo.commons_file}")
                    data = _download(client, photo.commons_url(IMAGE_WIDTH))
                    ext = Path(photo.commons_file).suffix.lower()
                    photo.image.save(f"{photo.organism_id}{ext}", ContentFile(data), save=False)
                photo.thumbnail.save(
                    f"{photo.organism_id}.jpg", _make_thumbnail(photo.image), save=False
                )
                photo.save()
        self.stdout.write(self.style.SUCCESS("Done"))


def _download(client: httpx.Client, url: str) -> bytes:
    """Fetch `url`, waiting out Wikimedia's 429s as its robot policy asks."""
    for _ in range(10):
        response = client.get(url)
        if response.status_code != 429:
            response.raise_for_status()
            return response.content
        time.sleep(int(response.headers.get("Retry-After", 10)))
    raise RuntimeError(f"Wikimedia kept rate-limiting {url}")


def _make_thumbnail(image) -> ContentFile:
    image.open("rb")
    with Image.open(image) as im:
        im = im.convert("RGB")
        im.thumbnail((THUMB_WIDTH, THUMB_WIDTH * 10))
        buf = BytesIO()
        im.save(buf, "JPEG", quality=85, optimize=True)
    image.close()
    return ContentFile(buf.getvalue())
