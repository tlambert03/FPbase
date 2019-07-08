from django.conf.urls import url
from favit import views

app_name = "favit"

urlpatterns = [
    url(r"^add-or-remove$", views.add_or_remove, name="add_or_remove"),
    url(r"^remove$", views.remove, name="remove"),
]
