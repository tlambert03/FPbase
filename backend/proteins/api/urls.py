from django.urls import path
from django.views.generic import TemplateView
from django.views.generic.base import RedirectView

from . import views

app_name = "api"

urlpatterns = [
    path("", TemplateView.as_view(template_name="pages/api.html"), name="api"),
    # /proteins/api
    path("spectrum/", views.SpectrumList.as_view(), name="get-spectrum"),
    path("spectrum/<int:pk>/", views.SpectrumDetail.as_view()),
    path("proteins/", views.ProteinListAPIView.as_view(), name="protein-api"),
    # path("proteins2/", views.ProteinListAPIView2.as_view(), name="protein-api2"),
    path(
        "proteins/table-data/",
        views.ProteinTableAPIView.as_view(),
        name="protein-table-api",
    ),
    path(
        "proteins/spectra/",
        views.ProteinSpectraListAPIView.as_view(),
        name="spectra-api",
    ),
    path(
        "proteins/basic/",
        views.BasicProteinListAPIView.as_view(),
        name="basic-protein-api",
    ),
    path("proteins/states/", views.StatesListAPIView.as_view(), name="states-api"),
    # /proteins/:slug/
    # re_path(
    #     r"^(?P<slug>[-\w]+)/$",
    #     views.ProteinRetrieveAPIView.as_view(),
    #     name="protein-api",
    # ),
    # non-normal endpoints
    path("proteins/spectraslugs/", RedirectView.as_view(url="/api/spectra-list/", permanent=True)),
    path("spectra-list/", views.spectra_list, name="spectra-list"),
    path("proteins/ocinfo/", RedirectView.as_view(url="/api/optical-configs-list/", permanent=True)),
    path("optical-configs-list/", views.optical_configs_list, name="ocinfo"),
]
