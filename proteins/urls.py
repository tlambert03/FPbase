from django.conf.urls import url
from django.views.generic import TemplateView

from . import views
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
                       message="You must be logged in to update protein information"),
        name='update'),
    url(r'^table/', views.protein_table, name='table'),
    url(r'^spectra/submit',
        login_required(views.SpectrumCreateView.as_view(),
                       message="You must be logged in to submit a new spectrum"),
        name='submit-spectra'),
    url(r'^spectra/(?P<slug>[-\w]+)', views.protein_spectra, name='spectra'),
    url(r'^spectra/', views.protein_spectra, name='spectra'),
    url(r'^spectra_csv/', views.spectra_csv, name='spectra_csv'),

    url(r'^chart/', TemplateView.as_view(template_name='ichart.html'), name='ichart'),
    url(r'^collections/(?P<owner>[-\w.]+)/?$', views.CollectionList.as_view(), name='collections'),
    url(r'^collections/', views.CollectionList.as_view(), name='collections'),
    url(r'^collection/(?P<pk>\d+)/$', views.CollectionDetail.as_view(), name='collection-detail'),
    url(r'^collection/(?P<pk>\d+)/update/',
        login_required(views.CollectionUpdateView.as_view(),
                       message="You must be logged in to update collections"),
        name='updatecollection'),
    url(r'^collection/(?P<pk>\d+)/delete/',
        login_required(views.CollectionDeleteView.as_view(),
                       message="You must be logged in to delete collections"),
        name='deletecollection'),
    url(r'^collection/create/',
        login_required(views.CollectionCreateView.as_view(),
                       message="You must be logged in to create a new collection"),
        name='newcollection'),

    url(r'^organism/(?P<pk>\d+)/$', views.OrganismDetailView.as_view(), name='organism-detail'),

    url(r'^protein/(?P<slug>[-\w]+)/$', views.ProteinDetailView.as_view(), name='protein-detail'),
    url(r'^protein/(?P<slug>[-\w]+)/bleach/$', views.protein_bleach_formsets, name='protein-bleach-form'),
    url(r'^protein/(?P<slug>[-\w]+)/rev/(?P<rev>\d+)$', views.ProteinDetailView.as_view(), name='protein-detail'),
    url(r'^protein/(?P<slug>[-\w]+)/ver/(?P<ver>\d+)$', views.ProteinDetailView.as_view(), name='protein-detail'),

    url(r'^protein-autocomplete/$', views.ProteinAutocomplete.as_view(), name='protein-autocomplete',),
    url(r'^state-autocomplete/$', views.StateAutocomplete.as_view(), name='state-autocomplete',),

    # /proteins/api
    url(r'^api/proteins/$', apiviews.ProteinListAPIView.as_view(), name='protein-api'),
    url(r'^api/proteins/spectraslugs/$', apiviews.spectraslugs, name='spectra-slugs'),
    url(r'^api/proteins/spectra/$', apiviews.ProteinSpectraListAPIView.as_view(), name='spectra-api'),
    url(r'^api/proteins/basic/$', apiviews.BasicProteinListAPIView.as_view(), name='basic-protein-api'),
    url(r'^api/proteins/states/$', apiviews.StatesListAPIView.as_view(), name='states-api'),
    # /proteins/api/:slug/
    url(r'^api/(?P<slug>[-\w]+)/$', apiviews.ProteinRetrieveUpdateDestroyAPIView.as_view(), name='protein-api'),

    # AJAX
    url(r'^ajax/add_taxonomy/$', views.add_organism, name='add_taxonomy'),
    url(r'^ajax/filter_import/(?P<brand>[-\w]+)$', views.filter_import, name='filter_import'),
    url(r'^ajax/add_protein_reference/(?P<slug>[-\w]+)/$', views.add_reference, name='add_protein_reference'),
    url(r'^ajax/admin_approve_protein/(?P<slug>[-\w]+)/$', views.approve_protein, name='admin_approve_protein'),
    url(r'^ajax/admin_revert_version/(?P<ver>\d+)$', views.revert_version, name='admin_revert_version'),
    url(r'^ajax/update_transitions/(?P<slug>[-\w]+)/$', views.update_transitions, name='update_transitions'),
    url(r'^ajax/validate_proteinname/$', views.validate_proteinname, name='validate_proteinname'),
    url(r'^ajax/validate_spectrumownername/$', views.similar_spectrum_owners, name='validate_spectrumownername'),
    url(r'^ajax/remove_from_collection/$', views.collection_remove, name='collection-remove'),
    url(r'^ajax/add_to_collection/$', views.add_to_collection, name='add_to_collection'),

]
