from django.db import models
from django.contrib.auth import get_user_model
from django.core.validators import URLValidator
User = get_user_model()


class Authorable(models.Model):
    created_by = models.ForeignKey(User, blank=True, null=True, related_name='%(class)s_author', on_delete=models.SET_NULL)
    updated_by = models.ForeignKey(User, blank=True, null=True, related_name='%(class)s_modifier', on_delete=models.SET_NULL)

    class Meta:
        abstract = True


class Product(models.Model):
    PRODUCT_LINKS = {
        'chroma': 'https://www.chroma.com/products/parts/*',
        'semrock': 'https://www.semrock.com/FilterDetails.aspx?id=*',
        'lumencor': 'http://lumencor.com/products/filters-for-spectra-x-light-engines/',
    }

    manufacturer = models.CharField(max_length=128, blank=True)
    part = models.CharField(max_length=128, blank=True)

    class Meta:
        abstract = True

    def get_absolute_url(self):
        part = self.part
        if self.manufacturer.lower() == 'chroma':
            part = part.replace('/', '-')

        try:
            url = self.PRODUCT_LINKS[self.manufacturer.lower()].replace('*', part)
            urlv = URLValidator()
            urlv(url)
            return url
        except Exception:
            return None
