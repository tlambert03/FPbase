from django.test import TestCase

from ..models import Protein, Spectrum, State


class TestProteinModel(TestCase):
    @classmethod
    def setUpTestData(cls):
        cls.protA: Protein = Protein.objects.create(name="ProteinA")
        cls.mprotB: Protein = Protein.objects.create(name="mProteinB")

        state = State.objects.create(
            ex_max=580,
            em_max=610,
            qy=0.22,
            pka=4.5,
            lifetime=1.4,
            protein=cls.protA,
        )

        Spectrum.objects.create(
            category=Spectrum.PROTEIN,
            subtype=Spectrum.EX,
            owner_state=state,
            data=[[431.0, 0.0039], [432.0, 0.0038], [433.0, 0.0042]],
        )
        Spectrum.objects.create(
            category=Spectrum.PROTEIN,
            subtype=Spectrum.EM,
            owner_state=state,
            data=[[550.0, 0.0012], [551.0, 0.0014], [552.0, 0.0015]],
        )

        cls.protA.default_state = state
        cls.protA.save()

    def test_protein_str(self):
        self.assertEqual(str(self.protA), "ProteinA")

    def test_get_absolute_url(self):
        # This will also fail if the urlconf is not defined.
        self.assertEqual(self.protA.get_absolute_url(), f"/protein/{self.protA.slug}/")

    def test_mless(self):
        # This will also fail if the urlconf is not defined.
        self.assertEqual(self.mprotB.mless, "ProteinB")

    def test_name(self):
        # This will also fail if the urlconf is not defined.
        self.assertEqual(str(self.protA), "ProteinA")

    def test_spectra_img(self):
        for ext in ("png", "svg", "tif", "pdf", "jpeg"):
            assert self.protA.spectra_img(fmt=ext)

    def test_spectrum_y_cache_invalidation(self):
        """Ensure cached y property is invalidated when change_y is called."""
        spectrum = Spectrum.objects.first()
        original_y = spectrum.y.copy()
        new_y = [0.5, 0.5, 0.5]
        spectrum.change_y(new_y)
        # The cached property should be invalidated and return the new values
        assert spectrum.y == new_y, f"Expected {new_y}, got {spectrum.y}"
        assert spectrum.y != original_y

    def test_spectrum_x_cache_invalidation(self):
        """Ensure cached x property is invalidated when change_x is called."""
        spectrum = Spectrum.objects.first()
        original_x = spectrum.x.copy()
        new_x = [500.0, 501.0, 502.0]
        spectrum.change_x(new_x)
        # The cached property should be invalidated and return the new values
        assert spectrum.x == new_x, f"Expected {new_x}, got {spectrum.x}"
        assert spectrum.x != original_x
