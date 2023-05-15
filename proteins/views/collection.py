import contextlib
import json

from django import forms
from django.contrib.auth.decorators import login_required
from django.core.exceptions import PermissionDenied
from django.core.mail import mail_admins
from django.http import (
    HttpResponseBadRequest,
    HttpResponseNotAllowed,
    HttpResponseRedirect,
    JsonResponse,
)
from django.shortcuts import get_object_or_404, render
from django.urls import resolve, reverse_lazy
from django.utils.text import slugify
from django.views.generic import (
    CreateView,
    DeleteView,
    DetailView,
    ListView,
    UpdateView,
)

from fpbase.util import is_ajax

from ..forms import CollectionForm
from ..models import ProteinCollection
from .mixins import OwnableObject


def serialized_proteins_response(queryset, format="json", filename="FPbase_proteins"):
    from proteins.api.serializers import ProteinSerializer

    ProteinSerializer.Meta.on_demand_fields = ()
    serializer = ProteinSerializer(queryset, many=True)
    if format == "csv":
        from django.http import StreamingHttpResponse
        from rest_framework_csv.renderers import CSVStreamingRenderer

        response = StreamingHttpResponse(CSVStreamingRenderer().render(serializer.data), content_type="text/csv")
        response["Content-Disposition"] = f'attachment; filename="{filename}.csv"'
    elif format == "json":
        response = JsonResponse(serializer.data, safe=False)
    return response


class CollectionList(ListView):
    def get_queryset(self):
        # get all collections for current user and all other non-private collections
        qs = ProteinCollection.objects.exclude(private=True).prefetch_related("owner")
        if self.request.user.is_authenticated:
            qs = qs | ProteinCollection.objects.filter(owner=self.request.user)
        if "owner" in self.kwargs:
            qs = qs.filter(owner__username=self.kwargs["owner"])
        return qs.order_by("-created")

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        if "owner" in self.kwargs:
            context["owner"] = self.kwargs["owner"]
        return context


class CollectionDetail(DetailView):
    queryset = ProteinCollection.objects.all().prefetch_related(
        "proteins",
        "proteins__states",
        "proteins__states__spectra",
        "proteins__default_state",
    )

    def get(self, request, *args, **kwargs):
        fmt = request.GET.get("format", "").lower()
        if fmt in ("json", "csv"):
            col = self.get_object()
            return serialized_proteins_response(col.proteins.all(), fmt, filename=slugify(col.name))
        return super().get(request, *args, **kwargs)

    def get_context_data(self, **kwargs):
        # Call the base implementation first to get a context
        context = super().get_context_data(**kwargs)
        context["isowner"] = self.request.user == self.object.owner

        _ids = []
        for prot in self.object.proteins.all():
            for state in prot.states.all():
                _ids.extend(sp.id for sp in state.spectra.all())

        context["spectra_ids"] = ",".join([str(i) for i in _ids])
        return context

    def render_to_response(self, *args, **kwargs):
        if not self.request.user.is_superuser and self.object.private and (self.object.owner != self.request.user):
            return render(self.request, "proteins/private_collection.html", {"foo": "bar"})
        return super().render_to_response(*args, **kwargs)


@login_required
def collection_remove(request):
    if not is_ajax(request):
        return HttpResponseNotAllowed([])
    try:
        protein = int(request.POST["target_protein"])
        collection = int(request.POST["target_collection"])
    except (KeyError, ValueError):
        return HttpResponseBadRequest()

    col = get_object_or_404(ProteinCollection, id=collection)

    if col.owner != request.user:
        return HttpResponseNotAllowed([])
    col.proteins.remove(protein)
    return JsonResponse({"status": "deleted"})


@login_required
def add_to_collection(request):
    if not is_ajax(request):
        return HttpResponseNotAllowed([])

    if request.method == "GET":
        qs = ProteinCollection.objects.filter(owner=request.user)
        widget = forms.Select(attrs={"class": "form-control custom-select", "id": "collectionSelect"})
        choicefield = forms.ChoiceField(choices=qs.values_list("id", "name"), widget=widget)

        members = []
        if request.GET.get("id"):
            with contextlib.suppress(Exception):
                qs = qs.filter(proteins=int(request.GET.get("id")))
                members = [(item.name, item.get_absolute_url()) for item in qs]
        response = {
            "widget": choicefield.widget.render("collectionChoice", ""),
            "members": json.dumps(members),
        }
        return JsonResponse(response)

    elif request.method == "POST":
        try:
            collection = ProteinCollection.objects.get(id=request.POST.get("collectionChoice"))
            collection.proteins.add(int(request.POST.get("protein")))
            status = "success"
        except Exception:
            status = "error"
        return JsonResponse({"status": status})

    return HttpResponseNotAllowed([])


class CollectionCreateView(OwnableObject, CreateView):
    model = ProteinCollection
    form_class = CollectionForm

    def get_form_kwargs(self):
        # add current user to the kwargs for the CollectionForm
        kwargs = super().get_form_kwargs()
        # the protein_search.html template has a form that can send a list of
        # proteins to create a new collection
        if len(self.request.POST.getlist("protein")):
            kwargs["proteins"] = self.request.POST.getlist("protein")
        # alternatively, a list of proteins from an existing collection can
        # be used to make a new collection
        elif self.request.POST.get("dupcollection", False):
            id = self.request.POST.get("dupcollection")
            kwargs["proteins"] = [p.id for p in ProteinCollection.objects.get(id=id).proteins.all()]
        return kwargs

    def form_valid(self, form):
        self.attach_owner(form)
        if getattr(form, "proteins", None):
            self.object.proteins.add(*form.proteins)
        if not self.request.user.is_staff:
            mail_admins(
                "Collection Created",
                f"User: {self.request.user.username}\nCollection: {self.object}\n"
                f"{self.request.build_absolute_uri(self.object.get_absolute_url())}",
                fail_silently=True,
            )
        return HttpResponseRedirect(self.get_success_url())

    def get_success_url(self):
        redirect_url = self.request.POST.get("next") or self.request.GET.get("next", None)
        try:
            # check that this is an internal redirection
            resolve(redirect_url)
        except Exception:
            redirect_url = None
        return redirect_url or super().get_success_url()

    def get_context_data(self, **kwargs):
        data = super().get_context_data(**kwargs)
        if self.request.POST.get("colname", False):
            data["colname"] = self.request.POST.get("colname")
        return data


class CollectionUpdateView(OwnableObject, UpdateView):
    model = ProteinCollection
    form_class = CollectionForm

    def dispatch(self, request, *args, **kwargs):
        if not self.get_object().has_change_permission(self.request):
            raise PermissionDenied
        return super().dispatch(request, *args, **kwargs)


class CollectionDeleteView(DeleteView):
    model = ProteinCollection
    success_url = reverse_lazy("proteins:collections")

    def dispatch(self, request, *args, **kwargs):
        if self.get_object().owner != self.request.user:
            raise PermissionDenied
        return super().dispatch(request, *args, **kwargs)

    def get_success_url(self):
        redirect_url = reverse_lazy("proteins:collections", kwargs={"owner": self.request.user})
        try:
            # check that this is an internal redirection
            resolve(redirect_url)
        except Exception:
            redirect_url = None
        return redirect_url or super().get_success_url()
