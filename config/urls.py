from django.conf import settings
from django.conf.urls.static import static
from django.contrib import admin
from django.contrib.sitemaps.views import sitemap
from django.urls import include, path, re_path
from django.views import defaults as default_views
from django.views.decorators.cache import cache_page
from django.views.decorators.csrf import csrf_exempt
from django.views.generic import TemplateView
from django.views.generic.base import RedirectView
from drf_spectacular.views import SpectacularAPIView, SpectacularSwaggerView
from graphene_django.views import GraphQLView

import fpbase.views
from fpbase.sitemaps import (
    AuthorsSitemap,
    MicroscopeSitemap,
    OrganismsSitemap,
    ProteinCollectionSitemap,
    ProteinSitemap,
    ReferencesSitemap,
    StaticSitemap,
)
from references.views import ReferenceListView

sitemaps = {
    "static": StaticSitemap(),
    "authors": AuthorsSitemap(),
    "references": ReferencesSitemap(),
    "proteins": ProteinSitemap(),
    "organisms": OrganismsSitemap(),
    "microscopes": MicroscopeSitemap(),
    "collections": ProteinCollectionSitemap(),
}

urlpatterns = [  # noqa: RUF005
    path("", fpbase.views.HomeView.as_view(), name="home"),
    path(
        "about/",
        TemplateView.as_view(template_name="pages/about.html"),
        name="about",
    ),
    # path('faq/', TemplateView.as_view(template_name='pages/faq.html'), name='faq'),
    path(
        "help/",
        RedirectView.as_view(url="https://help.fpbase.org", query_string=True),
        name="help",
    ),
    path("cite/", TemplateView.as_view(template_name="pages/cite.html"), name="cite"),
    path(
        "terms/",
        TemplateView.as_view(template_name="pages/terms.html"),
        name="terms",
    ),
    path(
        "privacy/",
        TemplateView.as_view(template_name="pages/terms.html"),
        name="privacy",
    ),
    path(
        "contributing/",
        TemplateView.as_view(template_name="pages/contributing.html"),
        name="contributing",
    ),
    path(
        "schema/",
        RedirectView.as_view(url="https://help.fpbase.org/schema/schema"),
        name="schema",
    ),
    path(
        "bleaching/",
        TemplateView.as_view(template_name="pages/bleaching.html"),
        name="bleaching",
    ),
    # path('mutations/', TemplateView.as_view(template_name='pages/mutations.html'), name='mutations'),
    re_path(
        r"^robots\.txt$",
        TemplateView.as_view(template_name="robots.txt", content_type="text/plain"),
        name="robots",
    ),
    re_path(
        r"^sitemap\.xml$",
        sitemap,
        {"sitemaps": sitemaps},
        name="django.contrib.sitemaps.views.sitemap",
    ),
    re_path(
        r"^googleaecf5301782589e7\.html$",
        TemplateView.as_view(template_name="googleaecf5301782589e7.html"),
        name="verification",
    ),
    path("contact/", fpbase.views.ContactView.as_view(), name="contact"),
    path(
        "thanks/",
        TemplateView.as_view(template_name="pages/thanks.html"),
        name="thanks",
    ),
    path("beta/", TemplateView.as_view(template_name="pages/beta.html"), name="beta"),
    # Django Admin, use {% url 'admin:index' %}
    re_path(settings.ADMIN_URL, admin.site.urls),
    # User management
    path("users/", include("fpbase.users.urls", namespace="users")),
    path("accounts/", include("allauth.urls")),
    # API base url
    # path("api/", include("proteins.api.urls", namespace="api")),
    # path('api/', TemplateView.as_view(template_name='pages/api.html'), name='api'),
    path("api/", include("config.api_router")),
    # api-auth for DRF
    path("api-auth/", include("rest_framework.urls", namespace="rest_framework")),
    # re_path(r"^api-docs/", include_docs_urls(title="FPbase API docs")),
    path("api/schema/", SpectacularAPIView.as_view(), name="api-schema"),
    path(
        "api/docs/",
        SpectacularSwaggerView.as_view(url_name="api-schema"),
        name="api-docs",
    ),
    # custom apps
    path("", include("proteins.urls")),  # NOTE: without $
    path("reference/", include("references.urls")),  # NOTE: without $
    path(
        "references/",
        cache_page(60 * 30)(ReferenceListView.as_view()),
        name="reference-list",
    ),
    path("fav/", include("favit.urls")),
    path("avatar/", include("avatar.urls")),
    re_path(r"^test500/", fpbase.views.test500),
    path("graphql/", csrf_exempt(GraphQLView.as_view(graphiql=True))),
    path("graphql/batch/", csrf_exempt(GraphQLView.as_view(batch=True))),
] + static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)

handler500 = "fpbase.views.server_error"

if settings.DEBUG:
    # This allows the error pages to be debugged during development, just visit
    # these url in browser to see how these error pages look like.

    urlpatterns += [
        path(
            "400/",
            default_views.bad_request,
            kwargs={"exception": Exception("Bad Request!")},
        ),
        path(
            "403/",
            default_views.permission_denied,
            kwargs={"exception": Exception("Permission Denied")},
        ),
        path(
            "404/",
            default_views.page_not_found,
            kwargs={"exception": Exception("Page not Found")},
        ),
        path("500/", fpbase.views.server_error),
        path("test/", fpbase.views.testview),
        path(
            "autocomplete/",
            TemplateView.as_view(template_name="pages/test_autocomplete.html"),
        ),
    ]
    if "debug_toolbar" in settings.INSTALLED_APPS:
        import debug_toolbar

        urlpatterns = [path("__debug__/", include(debug_toolbar.urls)), *urlpatterns]
