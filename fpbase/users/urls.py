from django.conf.urls import url
from django.contrib.admin.views.decorators import staff_member_required

from . import views

urlpatterns = [
    url(
        regex=r'^$',
        view=staff_member_required(views.UserListView.as_view()),
        name='list'
    ),
    url(
        regex=r'^~redirect/$',
        view=views.UserRedirectView.as_view(),
        name='redirect'
    ),
    url(
        regex=r'^(?P<username>[\w.@+-]+)/$',
        view=views.UserDetailView.as_view(),
        name='detail'
    ),
    url(
        regex=r'^~update/$',
        view=views.UserUpdateView.as_view(),
        name='update'
    ),
]
