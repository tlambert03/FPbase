from django.urls import path, re_path

from .views import (
    AuthorDetailView,
    ReferenceAutocomplete,
    ReferenceDetailView,
    add_excerpt,
)

app_name = "references"

urlpatterns = [
    path(
        "autocomplete/",
        ReferenceAutocomplete.as_view(create_field="doi"),
        name="reference-autocomplete",
    ),
    re_path(r"^author/(?P<pk>[-\w]+)/$", AuthorDetailView.as_view(), name="author-detail"),
    re_path(r"^ajax/add_excerpt/(?P<pk>[-\w]+)$", add_excerpt, name="add_excerpt"),
    re_path(r"^(?P<pk>[-\w\/\.]+)/$", ReferenceDetailView.as_view(), name="reference-detail"),
]
