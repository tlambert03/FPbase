from django.urls import path

from favit import views

app_name = "favit"

urlpatterns = [
    path("add-or-remove", views.add_or_remove, name="add_or_remove"),
    path("remove", views.remove, name="remove"),
]
