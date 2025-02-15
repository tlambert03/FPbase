import contextlib
import json
from textwrap import dedent

# from django.views.decorators.cache import cache_page
# from django.views.decorators.vary import vary_on_cookie
from django import forms
from django.conf import settings
from django.core.mail import EmailMessage
from django.http import Http404, HttpResponse, JsonResponse
from django.shortcuts import render
from django.template.defaultfilters import slugify
from django.urls import reverse_lazy
from django.views.decorators.clickjacking import xframe_options_exempt
from django.views.generic import CreateView

from fpbase.util import is_ajax, uncache_protein_page

from ..forms import SpectrumForm
from ..models import Filter, Protein, Spectrum, State
from ..util.importers import add_filter_to_database
from ..util.spectra import spectra2csv


# @cache_page(60 * 10)
# @vary_on_cookie
def protein_spectra(request, slug=None):
    """renders html for protein spectra page"""
    template = "spectra.html"

    if is_ajax(request) and slug is not None:
        """slug represents the slug of the OWNER of the spectra...
        for instance, a protein state.
        function returns an array containing ALL of the spectra
        belonging to that owner
        """
        owner = Spectrum.objects.filter_owner(slug)
        if owner.count():
            return JsonResponse({"spectra": [s.d3dict() for s in owner]})
        raise Http404

    return render(request, template)


@xframe_options_exempt
def protein_spectra_graph(request, slug=None):
    """Simple variant of protein_spectra only returns the graph.

    Used for embedding in iframes.
    """
    return render(request, "spectra_graph.html")


def spectrum_submitted(request):
    context = {"spectrum_name": request.session.get("spectrum_name", "")}
    return render(request, "spectrum_submitted.html", context)


class SpectrumCreateView(CreateView):
    model = Spectrum
    form_class = SpectrumForm

    def get_success_url(self, **kwargs):
        return reverse_lazy("proteins:spectrum_submitted")

    def get_initial(self):
        init = super().get_initial()
        if self.kwargs.get("slug", False):
            with contextlib.suppress(Exception):
                self.protein = Protein.objects.get(slug=self.kwargs.get("slug"))
                init["owner_state"] = self.protein.default_state
                init["category"] = Spectrum.PROTEIN

        return init

    def get_form_kwargs(self):
        kwargs = super().get_form_kwargs()
        kwargs["user"] = self.request.user
        return kwargs

    def get_form(self, form_class=None):
        form = super().get_form()

        if self.kwargs.get("slug", False):
            with contextlib.suppress(Exception):
                form.fields["owner_state"] = forms.ModelChoiceField(
                    required=True,
                    label="Protein (state)",
                    empty_label=None,
                    queryset=State.objects.filter(protein=self.protein).select_related("protein"),
                )
                form.fields["category"].disabled = True
        return form

    def form_valid(self, form):
        # Set the status to "pending" before saving the form
        if not self.request.user.is_staff:
            form.instance.status = Spectrum.STATUS.pending
        # This method is called when valid form data has been POSTed.
        # It should return an HttpResponse.
        response = super().form_valid(form)
        with contextlib.suppress(Exception):
            uncache_protein_page(self.object.owner_state.protein.slug, self.request)

        if not self.request.user.is_staff:
            body = f"""
            A new spectrum has been submitted to FPbase.

            Admin URL: {self.request.build_absolute_uri(form.instance.get_admin_url())}
            User: {self.request.user}
            Data:

            {form.cleaned_data}
            """
            EmailMessage(
                subject=f"[FPbase] Spectrum needs validation: {form.cleaned_data['owner']}",
                body=dedent(body),
                to=[a[1] for a in settings.ADMINS],
                headers={"X-Mailgun-Track": "no"},
            ).send()
        self.request.session["spectrum_name"] = self.object.name
        return response

    def get_context_data(self, **kwargs):
        # Call the base implementation first to get a context
        context = super().get_context_data(**kwargs)
        if hasattr(self, "protein"):
            context["protein"] = self.protein
        return context


def spectra_csv(request):
    try:
        idlist = [int(x) for x in request.GET.get("q", "").split(",") if x]
        spectralist = Spectrum.objects.filter(id__in=idlist)
        if spectralist:
            return spectra2csv(spectralist)
    except Exception:
        return HttpResponse("malformed spectra csv request")
    else:
        return HttpResponse("malformed spectra csv request")


def filter_import(request, brand):
    part = request.POST["part"]
    new_objects = []
    errors = []
    response = {"status": 0}

    with contextlib.suppress(Filter.DoesNotExist):
        Filter.objects.get(slug=slugify(f"{brand} {part}"))
        response["message"] = f"{part} is already in the database"
        return JsonResponse(response)

    try:
        new_objects, errors = add_filter_to_database(brand, part, request.user)
    except Exception as e:
        response["message"] = str(e)

    if new_objects:
        spectrum = new_objects[0]
        response = {
            "status": 1,
            "objects": spectrum.name,
            "spectra_options": json.dumps(
                {
                    "category": spectrum.category,
                    "subtype": spectrum.subtype,
                    "slug": spectrum.owner.slug,
                    "name": spectrum.owner.name,
                }
            ),
        }
    elif errors:
        with contextlib.suppress(Exception):
            if errors[0][1].as_data()["owner"][0].code == "owner_exists":
                response["message"] = f"{part} already appears to be imported"
    return JsonResponse(response)
