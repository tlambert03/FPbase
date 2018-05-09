from django.conf import settings
from django.conf.urls import include, url
from django.conf.urls.static import static
from django.contrib import admin
from django.views.generic import TemplateView
from django.views import defaults as default_views
from django.contrib.sitemaps.views import sitemap
from rest_framework.documentation import include_docs_urls
from fpbase.sitemaps import ProteinSitemap, OrganismsSitemap, StaticSitemap
import fpbase.views

sitemaps = {
    'static': StaticSitemap(),
    'proteins': ProteinSitemap(),
    # 'protstat': ProteinStaticSitemap(),
    'organisms': OrganismsSitemap(),
}


urlpatterns = [
    url(r'^$', TemplateView.as_view(template_name='pages/home.html'), name='home'),
    url(r'^about/$', TemplateView.as_view(template_name='pages/about.html'), name='about'),
    url(r'^terms/$', TemplateView.as_view(template_name='pages/terms.html'), name='terms'),
    url(r'^privacy/$', TemplateView.as_view(template_name='pages/terms.html'), name='privacy'),
    url(r'^contributing/$', TemplateView.as_view(template_name='pages/contributing.html'), name='contributing'),
    url(r'^schema/$', TemplateView.as_view(template_name='pages/schema.html'), name='schema'),
    url(r'^robots\.txt$', TemplateView.as_view(template_name="robots.txt", content_type="text/plain"), name="robots"),
    url(r'^sitemap\.xml$', sitemap, {'sitemaps': sitemaps}, name='django.contrib.sitemaps.views.sitemap'),
    url(r'^googleaecf5301782589e7\.html$', TemplateView.as_view(template_name="googleaecf5301782589e7.html"), name="verification"),

    url(r'^contact/$', fpbase.views.ContactView.as_view(), name='contact'),
    url(r'^thanks/$', TemplateView.as_view(template_name='pages/thanks.html'), name='thanks'),
    url(r'^beta/$', TemplateView.as_view(template_name='pages/beta.html'), name='beta'),

    # Django Admin, use {% url 'admin:index' %}
    url(settings.ADMIN_URL, admin.site.urls),

    # User management
    url(r'^users/', include('fpbase.users.urls', namespace='users')),
    url(r'^accounts/', include('allauth.urls')),
    url(r'^api/$', TemplateView.as_view(template_name='pages/api.html'), name='api'),


    # Your stuff: custom urls includes go here

    # api-auth for DRF
    url(r'^api-auth/', include('rest_framework.urls', namespace='rest_framework')),
    url(r'^api-docs/', include_docs_urls(title='FPbase API docs')),

    # custom apps
    url(r'^', include('proteins.urls')),  # NOTE: without $
    url(r'^reference/', include('references.urls')),  # NOTE: without $

    url(r'^fav/', include('favit.urls')),

] + static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)


if settings.DEBUG:
    # This allows the error pages to be debugged during development, just visit
    # these url in browser to see how these error pages look like.
    urlpatterns += [
        url(r'^400/$', default_views.bad_request, kwargs={'exception': Exception('Bad Request!')}),
        url(r'^403/$', default_views.permission_denied, kwargs={'exception': Exception('Permission Denied')}),
        url(r'^404/$', default_views.page_not_found, kwargs={'exception': Exception('Page not Found')}),
        url(r'^500/$', default_views.server_error),
    ]
    if 'debug_toolbar' in settings.INSTALLED_APPS:
        import debug_toolbar
        urlpatterns = [
            url(r'^__debug__/', include(debug_toolbar.urls)),
        ] + urlpatterns
