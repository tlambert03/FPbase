from django.conf.urls import url
from django.views.generic import TemplateView

from . import views
from .views import ProteinDetailView
from .api import views as apiviews
from fpbase.decorators import login_required_message_and_redirect as login_required


app_name = 'proteins'

urlpatterns = [
    # detail view: /:slug
    url(r'^search/', views.protein_search, name='search'),
    url(r'^submit/',
        login_required(views.ProteinCreateView.as_view(),
                       message="You must be logged in to submit a new protein"),
        name='submit'),
    url(r'^protein/(?P<slug>[-\w]+)/update/',
        login_required(views.ProteinUpdateView.as_view(),
                       message="You must be logged in to submit a new protein"),
        name='update'),
    url(r'^table/', views.protein_table, name='table'),
    url(r'^chart/', TemplateView.as_view(template_name='ichart.html'), name='ichart'),

    url(r'^protein/(?P<slug>[-\w]+)/$', ProteinDetailView.as_view(), name='protein-detail'),

    # /proteins/api
    url(r'^api/$', apiviews.ProteinListCreateAPIView.as_view(), name='protein-api'),
    url(r'^api/basic/$', apiviews.BasicProteinListCreateAPIView.as_view(), name='basic-protein-api'),
    # /proteins/api/:slug/
    url(r'^api/(?P<slug>[-\w]+)/$', apiviews.ProteinRetrieveUpdateDestroyAPIView.as_view(), name='protein-api')

]