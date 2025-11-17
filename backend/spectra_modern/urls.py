"""URL configuration for modern spectrum submission."""

from django.urls import path

from . import views

app_name = "spectra_modern"

urlpatterns = [
    path("submit/", views.ModernSpectrumCreateView.as_view(), name="submit"),
    path("submitted/", views.spectrum_submitted, name="submitted"),
    path("preview/", views.spectrum_preview, name="preview"),
    path("api/state-search/", views.state_search, name="state_search"),
]
