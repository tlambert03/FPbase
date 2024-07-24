import proteins.api.urls
from django.conf import settings
from fpbase.users.api.views import UserViewSet
from rest_framework.routers import DefaultRouter, SimpleRouter

if settings.DEBUG:
    router = DefaultRouter()
else:
    router = SimpleRouter()

router.register("users", UserViewSet)

app_name = "api"


urlpatterns = proteins.api.urls.urlpatterns + router.urls
