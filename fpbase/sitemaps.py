from django.contrib.sitemaps import Sitemap
from django.urls import reverse
from proteins.models import Protein, Organism
from references.models import Author


class ProteinSitemap(Sitemap):
    changefreq = "weekly"
    priority = 0.7

    def items(self):
        return Protein.visible.all()

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


class StaticSitemap(Sitemap):
    priority = 0.8
    changefreq = 'weekly'

    # The below method returns all urls defined in urls.py file
    def items(self):
        from config.urls import urlpatterns as homeUrls
        return [url.name for url in homeUrls if hasattr(url, 'name') and url.name]

    def location(self, item):
        return reverse(item)
