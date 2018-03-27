from django.test import TestCase
from ..models import Protein


class TestProteinModel(TestCase):

    @classmethod
    def setUpTestData(cls):
        cls.protA = Protein.objects.create(name='ProteinA')
        cls.protB = Protein.objects.create(name='ProteinB')

    def test_protein_str(self):
        self.assertEquals(str(self.protA), 'ProteinA')

    def test_get_absolute_url(self):
        # This will also fail if the urlconf is not defined.
        self.assertEquals(self.protA.get_absolute_url(),
                          '/protein/{}/'.format(self.protA.slug))
