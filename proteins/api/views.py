from rest_framework.generics import (
    ListAPIView,
    ListCreateAPIView,
    RetrieveUpdateDestroyAPIView
)
from rest_framework.permissions import IsAuthenticated, IsAdminUser, AllowAny

from ..models import Protein
from .serializers import ProteinSerializer, BasicProteinSerializer

from django.utils.decorators import method_decorator
from django.views.decorators.cache import cache_page


class BasicProteinListCreateAPIView(ListAPIView):
    queryset = Protein.objects.filter(switch_type=Protein.BASIC).select_related('default_state')
    permission_classes = (AllowAny, )
    serializer_class = BasicProteinSerializer
    lookup_field = 'slug'  # Don't use Protein.id!

    @method_decorator(cache_page(60 * 15))
    def dispatch(self, *args, **kwargs):
        return super(BasicProteinListCreateAPIView, self).dispatch(*args, **kwargs)


class ProteinListCreateAPIView(ListCreateAPIView):
    queryset = Protein.objects.all().prefetch_related('states', 'transitions')
    permission_classes = (IsAuthenticated, )
    serializer_class = ProteinSerializer
    lookup_field = 'slug'  # Don't use Protein.id!


class ProteinRetrieveUpdateDestroyAPIView(RetrieveUpdateDestroyAPIView):
    queryset = Protein.objects.all()
    permission_classes = (IsAdminUser, )
    serializer_class = ProteinSerializer
    lookup_field = 'slug'  # Don't use Protein.id
