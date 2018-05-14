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
        'thermofisher': 'https://www.thermofisher.com/us/en/home/life-science/cell-analysis/fluorophores/*.html'
    }

    manufacturer = models.CharField(max_length=128, blank=True)
    part = models.CharField(max_length=128, blank=True)

    class Meta:
        abstract = True

    def get_absolute_url(self):
        part = self.part
        if not part:
            return None

        if self.manufacturer.lower() == 'chroma':
            part = part.replace('/', '-')

        try:
            if 'thermo' in self.manufacturer.lower():
                goodParts = ['bodipy-fl', 'coumarin', 'cy3-dye', 'cy5-dye', 'fluorescein',
                             'oregon-green', 'pacific-blue-dye', 'pacific-green-dye',
                             'pacific-orange-dye', 'tritc-dye', 'texas-red', 'dapi-stain',
                             'propidium-iodide', 'styo-9', 'sytox-green-stain', 'to-pro-3']
                alexas = (350, 405, 488, 532, 546, 555, 568, 594, 647, 680, 750)
                goodParts += ['alexa-fluor-%s' % a for a in alexas]
                goodParts += ['qdot-%s' % q for q in (525, 565, 605, 655, 705, 800)]
                assert part in goodParts, 'Thermo URL does not exist'
                url = self.PRODUCT_LINKS['thermofisher'].replace('*', part)
            else:
                url = self.PRODUCT_LINKS[self.manufacturer.lower()].replace('*', part)
            urlv = URLValidator()
            urlv(url)
            return url
        except Exception:
            return None
