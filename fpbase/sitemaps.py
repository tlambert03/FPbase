from django.contrib.sitemaps import Sitemap
from django.urls import reverse

from proteins.models import Microscope, Organism, Protein, ProteinCollection
from references.models import Author, Reference


class ProteinCollectionSitemap(Sitemap):
    changefreq = "weekly"
    priority = 0.7

    def items(self):
        return ProteinCollection.objects.exclude(private=True)

    def lastmod(self, item):
        return item.modified


class ProteinSitemap(Sitemap):
    changefreq = "weekly"
    priority = 0.7

    def items(self):
        return Protein.visible.all()

    def lastmod(self, item):
        return item.modified


class MicroscopeSitemap(Sitemap):
    changefreq = "weekly"
    priority = 0.7

    def items(self):
        return Microscope.objects.all()

    def lastmod(self, item):
        return item.modified


class OrganismsSitemap(Sitemap):
    changefreq = "weekly"
    priority = 0.6

    def items(self):
        return Organism.objects.all()

    def lastmod(self, item):
        return item.modified


class AuthorsSitemap(Sitemap):
    changefreq = "weekly"
    priority = 0.6

    def items(self):
        return Author.objects.all()

    def lastmod(self, item):
        return item.modified


class ReferencesSitemap(Sitemap):
    changefreq = "weekly"
    priority = 0.6

    def items(self):
        return Reference.objects.all()

    def lastmod(self, item):
        return item.modified


class StaticSitemap(Sitemap):
    priority = 0.8
    changefreq = "weekly"

    # The below method returns all urls defined in urls.py file
    def items(self):
        from config.urls import urlpatterns as homeUrls

        u = [url.name for url in homeUrls if hasattr(url, "name") and url.name]
        protUrls = [
            "search",
            "submit",
            "table",
            "submit-spectra",
            "spectra",
            "ichart",
            "collections",
            "microscopes",
            "fret",
            "lineage",
            "compare",
            "organisms",
            "activity",
            "blast",
        ]
        for url in protUrls:
            try:
                reverse("proteins:" + url)
                u.append("proteins:" + url)
            except Exception:
                pass
        return u

    def location(self, item):
        return reverse(item)
