from django.db import models
from django.contrib.auth import get_user_model
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
    url = models.URLField(blank=True)

    class Meta:
        abstract = True

    def save(self, *args, **kwargs):
        if self.part and (self.manufacturer.lower() in self.PRODUCT_LINKS):
            part = self.part
            if self.manufacturer.lower() == 'chroma':
                part = part.replace('/', '-')
            self.url = self.PRODUCT_LINKS[self.manufacturer.lower()].replace('*', part)
        super().save(*args, **kwargs)

    def get_absolute_url(self):
        if self.url:
            return self.url
