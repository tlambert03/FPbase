from django.conf.urls import url
from .views import ReferenceDetailView, AuthorDetailView

app_name = 'references'

urlpatterns = [
	url(r'^(?P<pk>[-\w]+)/$', ReferenceDetailView.as_view(), name='reference-detail'),
	url(r'^author/(?P<pk>[-\w]+)/$', AuthorDetailView.as_view(), name='author-detail'),
]