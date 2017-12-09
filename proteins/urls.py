from django.conf.urls import url
from django.views.generic import TemplateView

from . import views
from .views import ProteinDetailView
from .api import views as apiviews

app_name = 'proteins'

urlpatterns = [
	# detail view: /:slug
	url(r'^search/', views.search, name='search'),
	url(r'^submit/', views.submit),
	url(r'^table/', views.protein_table, name='table'),
	url(r'^chart/', TemplateView.as_view(template_name='ichart.html'), name='ichart'),

	url(r'^protein/(?P<slug>[-\w]+)/$', ProteinDetailView.as_view(), name='protein-detail'),

	# /proteins/api
	url(r'^api/$', apiviews.ProteinListCreateAPIView.as_view(), name='protein-api'),
	url(r'^api/basic/$', apiviews.BasicProteinListCreateAPIView.as_view(), name='basic-protein-api'),
	# /proteins/api/:slug/
	url(r'^api/(?P<slug>[-\w]+)/$', apiviews.ProteinRetrieveUpdateDestroyAPIView.as_view(), name='protein-api')

]