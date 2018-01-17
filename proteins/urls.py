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
    url(r'^chart/', TemplateView.as_view(template_name='ichart.html'), name='ichart'),
    url(r'^collections/(?P<owner>[-\w]+)/?$', views.CollectionList.as_view(), name='collections'),
    url(r'^collections/', views.CollectionList.as_view(), name='collections'),
    url(r'^collection/(?P<pk>\d+)/$', views.CollectionDetail.as_view(), name='collection-detail'),
    url(r'^collection/(?P<pk>\d+)/update/',
        login_required(views.CollectionUpdateView.as_view(),
        message="You must be logged in to delete collections"),
        name='updatecollection'),
    url(r'^newcollection/',
        login_required(views.CollectionCreateView.as_view(),
        message="You must be logged in to create a new collection"),
        name='newcollection'),


    url(r'^ajax/add_protein_reference/(?P<slug>[-\w]+)/$', views.add_reference, name='add_protein_reference'),
    url(r'^ajax/admin_approve_protein/(?P<slug>[-\w]+)/$', views.approve_protein, name='admin_approve_protein'),
    url(r'^ajax/admin_revert_version/(?P<ver>\d+)$', views.revert_version, name='admin_revert_version'),
    url(r'^ajax/update_transitions/(?P<slug>[-\w]+)/$', views.update_transitions, name='update_transitions'),
    url(r'^ajax/validate_proteinname/$', views.validate_proteinname, name='validate_proteinname'),
    url(r'^ajax/remove_from_collection/$', views.collection_remove, name='collection-remove'),
    url(r'^ajax/add_to_collection/$', views.add_to_collection, name='add_to_collection'),

    url(r'^protein/(?P<slug>[-\w]+)/$', views.ProteinDetailView.as_view(), name='protein-detail'),
    url(r'^protein/(?P<slug>[-\w]+)/rev/(?P<rev>\d+)$', views.ProteinDetailView.as_view(), name='protein-detail'),
    url(r'^protein/(?P<slug>[-\w]+)/ver/(?P<ver>\d+)$', views.ProteinDetailView.as_view(), name='protein-detail'),

    # /proteins/api
    url(r'^api/$', apiviews.ProteinListCreateAPIView.as_view(), name='protein-api'),
    url(r'^api/basic/$', apiviews.BasicProteinListCreateAPIView.as_view(), name='basic-protein-api'),
    # /proteins/api/:slug/
    url(r'^api/(?P<slug>[-\w]+)/$', apiviews.ProteinRetrieveUpdateDestroyAPIView.as_view(), name='protein-api')

]