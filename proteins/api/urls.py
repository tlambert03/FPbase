from django.urls import path
from django.views.generic import TemplateView

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
    path("proteins/spectraslugs/", views.spectraslugs, name="spectra-slugs"),
    path("proteins/ocinfo/", views.ocinfo, name="ocinfo"),
]
