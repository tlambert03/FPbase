from django.conf.urls import url
from django.urls import path
from django.views.generic import TemplateView
from . import views

app_name = 'api'

urlpatterns = [
    url(r'^$', TemplateView.as_view(template_name='pages/api.html'), name='api'),
    # /proteins/api
    url(r'^spectrum/$', views.SpectrumList.as_view(), name='get-spectrum'),
    path('spectrum/<int:pk>/', views.SpectrumDetail.as_view()),
    url(r'^proteins/$', views.ProteinListAPIView.as_view(), name='protein-api'),
    url(r'^proteins2/$', views.ProteinListAPIView2.as_view(), name='protein-api2'),
    url(r'^proteins/spectraslugs/$', views.spectraslugs, name='spectra-slugs'),
    url(r'^proteins/spectra/$', views.ProteinSpectraListAPIView.as_view(), name='spectra-api'),
    url(r'^proteins/basic/$', views.BasicProteinListAPIView.as_view(), name='basic-protein-api'),
    url(r'^proteins/states/$', views.StatesListAPIView.as_view(), name='states-api'),
    # /proteins/:slug/
    url(r'^(?P<slug>[-\w]+)/$', views.ProteinRetrieveAPIView.as_view(), name='protein-api'),

]

