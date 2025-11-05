import contextlib
import io
import json
import logging
import operator
from functools import reduce
from typing import TYPE_CHECKING, cast

import django.forms
import django.forms.formsets
import reversion
from django.apps import apps
from django.conf import settings
from django.contrib import messages
from django.contrib.admin.views.decorators import staff_member_required
from django.contrib.auth.decorators import login_required
from django.core.mail import mail_admins, mail_managers
from django.db import transaction
from django.db.models import Case, Count, F, Func, Prefetch, Q, Value, When
from django.forms.models import BaseInlineFormSet, modelformset_factory
from django.http import (
    Http404,
    HttpResponse,
    HttpResponseBadRequest,
    HttpResponseNotAllowed,
    HttpResponseRedirect,
    JsonResponse,
)
from django.shortcuts import get_object_or_404, redirect, render
from django.utils.decorators import method_decorator
from django.utils.html import strip_tags
from django.utils.safestring import mark_safe
from django.utils.text import slugify
from django.views.decorators.cache import cache_page
from django.views.decorators.vary import vary_on_cookie
from django.views.generic import CreateView, DetailView, ListView, UpdateView, base
from reversion.models import Revision, Version

from fpbase.util import is_ajax, uncache_protein_page
from proteins.extrest.entrez import get_cached_gbseqs
from proteins.extrest.ga import cached_ga_popular
from proteins.util.helpers import link_excerpts, most_favorited
from proteins.util.maintain import check_lineages, suggested_switch_type
from proteins.util.spectra import spectra2csv
from references.models import Reference  # breaks application modularity

from ..forms import (
    BleachComparisonForm,
    BleachMeasurementForm,
    LineageFormSet,
    ProteinForm,
    StateFormSet,
    StateTransitionFormSet,
    bleach_items_formset,
)
from ..models import BleachMeasurement, Excerpt, Organism, Protein, Spectrum, State

if TYPE_CHECKING:
    import maxminddb

    from proteins.forms.forms import BaseStateFormSet

logger = logging.getLogger(__name__)


def check_switch_type(object, request):
    suggested = suggested_switch_type(object)
    if suggested and object.switch_type != suggested:
        disp = dict(Protein.SWITCHING_CHOICES).get(suggested).lower()
        actual = object.get_switch_type_display().lower()
        msg = (
            "<i class='fa fa-exclamation-circle mr-2'></i><strong>Warning:</strong> "
            + "Based on the number of states and transitions currently assigned "
            + f"to this protein, it appears to be a {disp} protein; "
            + f"however, it has been assigned a switching type of {actual}. "
            + "Please confirm that the switch type, states, and transitions are correct."
        )
        messages.add_message(request, messages.WARNING, msg)


def form_changes(form, pre=""):
    if not form.has_changed():
        return []
    changes = []
    for field_name in form.changed_data:
        old = form.initial.get(field_name)
        new = form.cleaned_data.get(field_name)
        chg = f"{pre}{form.fields.get(field_name).label}: {old} -> {new}"
        changes.append(chg)
    return changes


def formset_changes(formset):
    if not formset.has_changed():
        return []

    model = formset.model.__name__
    changes = [f"Deleted {model} {obj}" for obj in formset.deleted_objects]
    changes.extend(f"Created {model} {obj}" for obj in formset.new_objects)
    for obj, _fields in formset.changed_objects:
        form = next(_form for _form in formset.forms if _form.instance == obj)
        frm_chg = form_changes(form)
        changes.append(f"Changed {model} {obj}: ({'; '.join(frm_chg)})")
    return changes


def get_form_changes(*forms):
    """accepts a list of forms or formsets and returns a list of all the fields
    that were added, deleted, or changed"""

    changes = []
    for form in forms:
        if isinstance(form, django.forms.BaseForm):
            changes.extend(form_changes(form))
        elif isinstance(form, django.forms.formsets.BaseFormSet):
            changes.extend(formset_changes(form))
    return changes


class _RollBackRevisionView(Exception):
    def __init__(self, response):
        self.response = response


def maxmind_db() -> str:
    """Create and return a temporary file containing the MaxMind database.

    Uses Django cache to store the database path across workers.
    The file persists on disk but the path is cached for 24 hours.
    """
    import io
    import os
    import tarfile
    import tempfile

    import requests
    from django.conf import settings
    from django.core.cache import cache

    # Try to get cached path first
    cache_key = "maxmind_db_path"
    cached_path = cache.get(cache_key)
    if cached_path and os.path.exists(cached_path):
        return cached_path

    # Download and cache the database
    url = "https://download.maxmind.com/app/geoip_download?edition_id=GeoLite2-Country&license_key={}&suffix=tar.gz"
    url = url.format(settings.MAXMIND_API_KEY)
    response = requests.get(url)
    response.raise_for_status()
    with tarfile.open(fileobj=io.BytesIO(response.content), mode="r:gz") as tar:
        for member in tar.getmembers():
            if member.name.endswith(".mmdb"):
                mmdb_file = tar.extractfile(member)
                if mmdb_file is not None:
                    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".mmdb")
                    tmp.write(mmdb_file.read())
                    tmp.close()
                    # Cache the path for 24 hours
                    cache.set(cache_key, tmp.name, 60 * 60 * 24)
                    return tmp.name
    return ""


def maxmind_reader() -> "maxminddb.Reader | None":
    """Get MaxMind database reader.

    Uses Django cache to minimize memory usage. The reader is cached
    for 1 hour to balance memory usage and performance.
    """
    from django.core.cache import cache
    from maxminddb import open_database

    cache_key = "maxmind_reader"
    reader = cache.get(cache_key)
    if reader is not None:
        return reader

    try:
        if db := maxmind_db():
            reader = open_database(db)
            # Cache reader for 1 hour
            cache.set(cache_key, reader, 60 * 60)
            return reader
    except Exception:
        pass
    return None


def get_country_code(request) -> str:
    # Definitely should be used inside a try/exc block
    if reader := maxmind_reader():
        x_forwarded_for = request.headers.get("x-forwarded-for")
        if x_forwarded_for:
            ip = x_forwarded_for.split(",")[0]
        else:
            ip = request.META.get("REMOTE_ADDR")
        if response := reader.get(ip):
            return str(response["country"]["iso_code"])  # pyright: ignore[reportIndexIssue]
    return ""


class ProteinDetailView(DetailView):
    """renders html for single protein page"""

    queryset = (
        Protein.objects.annotate(has_spectra=Count("states__spectra"))
        .prefetch_related(
            "states__spectra",  # Prefetch spectra to avoid N+1 queries
            "states",
            "excerpts__reference",
            "oser_measurements__reference",
        )
        .select_related("primary_reference")
    )

    def dispatch(self, *args, **kwargs):
        return super().dispatch(*args, **kwargs)

    # Only enable caching in production (when DEBUG=False)
    if not settings.DEBUG:
        dispatch = method_decorator(cache_page(60 * 30))(dispatch)
        dispatch = method_decorator(vary_on_cookie)(dispatch)

    def version_view(self, request, version, *args, **kwargs):
        try:
            with transaction.atomic(using=version.db):
                # Revert the revision.
                version.revision.revert(delete=True)
                # Run the normal changeform view.
                self.object = self.get_object()
                context = self.get_context_data(object=self.object)
                context["version"] = version
                response = self.render_to_response(context)
                response.render()  # eager rendering of response is necessary before db rollback
                raise _RollBackRevisionView(response)
        except _RollBackRevisionView as ex:
            return ex.response

    def get_object(self, queryset=None):
        if queryset is None:
            queryset = self.get_queryset()
        try:
            obj = queryset.get(slug=self.kwargs.get("slug"))
        except Protein.DoesNotExist:
            try:
                obj = queryset.get(uuid=self.kwargs.get("slug", "").upper())
            except Protein.DoesNotExist as e:
                raise Http404("No protein found matching this query") from e
        if obj.status == "hidden" and obj.created_by != self.request.user and not self.request.user.is_staff:
            raise Http404("No protein found matching this query")
        return obj

    def get(self, request, *args, **kwargs):
        if "rev" in kwargs:
            try:
                rev = int(kwargs["rev"])  # has to be int or indexing will fail
            except Exception:
                rev = 0
            if rev > 0:
                versions = Version.objects.get_for_object(self.get_object())
                version = versions[min(versions.count() - 1, rev)]
                return self.version_view(request, version, *args, **kwargs)
        elif "ver" in kwargs:
            version = get_object_or_404(Version, id=kwargs["ver"])
            if int(version.object_id) == self.get_object().id:
                return self.version_view(request, version, *args, **kwargs)
            # TODO:  ELSE WHAT??
        try:
            return super().get(request, *args, **kwargs)
        except Http404:
            from django.contrib.postgres.fields import ArrayField
            from django.db import models

            name = slugify(self.kwargs.get(self.slug_url_kwarg))
            aliases_lower = Func(Func(F("aliases"), function="unnest"), function="LOWER")
            remove_space = Func(aliases_lower, Value(" "), Value("-"), function="replace")
            final = Func(
                remove_space,
                Value("."),
                Value(""),
                function="replace",
                output_field=ArrayField(models.CharField(max_length=200)),
            )
            d = dict(Protein.objects.annotate(aka=final).values_list("aka", "id"))
            if name in d:
                obj = Protein.objects.get(id=d[name])
                messages.add_message(
                    self.request,
                    messages.INFO,
                    f"The URL {self.request.get_full_path()} was not found.  You have been forwarded here",
                )
                return HttpResponseRedirect(obj.get_absolute_url())
            raise

    def get_context_data(self, **kwargs):
        data = super().get_context_data(**kwargs)
        if self.object.status != "approved":
            data["last_approved"] = self.object.last_approved_version()

        similar = Protein.visible.filter(name__iexact=f"m{self.object.name}")
        similar = similar | Protein.visible.filter(name__iexact=f"monomeric{self.object.name}")
        similar = similar | Protein.visible.filter(name__iexact=self.object.name.lstrip("m"))
        similar = similar | Protein.visible.filter(name__iexact=self.object.name.lower().lstrip("td"))
        similar = similar | Protein.visible.filter(name__iexact=f"td{self.object.name}")
        data["similar"] = similar.exclude(id=self.object.id)
        spectra = [sp for state in self.object.states.all() for sp in state.spectra.all()]

        data["spectra_ids"] = ",".join([str(sp.id) for sp in spectra])
        data["hidden_spectra"] = ",".join([str(sp.id) for sp in spectra if sp.subtype in ("2p")])

        # put links in excerpts
        data["excerpts"] = link_excerpts(self.object.excerpts.all(), self.object.name, self.object.aliases)

        # Add country code to context
        try:
            data["country_code"] = get_country_code(self.request)
        except Exception:
            data["country_code"] = ""

        # Serialize PDB IDs as JSON for JavaScript
        if self.object.pdb:
            data["pdb_ids_json"] = json.dumps(self.object.pdb)

        return data


class ProteinCreateUpdateMixin:
    def get_form_type(self):
        return self.request.resolver_match.url_name

    def form_valid(self, form: ProteinForm):
        # This method is called when valid form data has been POSTed.
        context: dict = self.get_context_data()
        states = cast("BaseStateFormSet", context["states"])
        lineage = cast("BaseInlineFormSet", context["lineage"])

        if not (states.is_valid() and lineage.is_valid()):
            context |= {"states": states, "lineage": lineage}
            return self.render_to_response(context)

        # check if another protein already has the sequence that would be
        if not form.cleaned_data.get("seq"):
            for lform in lineage.forms:
                mut = lform.cleaned_data.get("mutation")
                par = lform.cleaned_data.get("parent")
                if par and mut:
                    seq = par.protein.seq.mutate(mut)
                    with contextlib.suppress(Protein.DoesNotExist):
                        prot = Protein.objects.get(seq__iexact=str(seq))
                        msg = mark_safe(
                            f'<a href="{prot.get_absolute_url()}" style="text-decoration: '
                            + f'underline;">{prot.name}</a> already has the sequence that'
                            + " would be generated by this parent & mutation"
                        )
                        lform.add_error(None, msg)
                        context |= {"states": states, "lineage": lineage}
                        return self.render_to_response(context)

        with transaction.atomic():
            with reversion.create_revision():
                self.object = form.save()
                self.object.primary_reference = (
                    Reference.objects.get_or_create(doi=doi.lower())[0]
                    if (doi := form.cleaned_data.get("reference_doi"))
                    else None
                )

                states.instance = self.object
                saved_states = states.save(commit=False)
                for s in saved_states:
                    if not s.created_by:
                        s.created_by = self.request.user
                    s.updated_by = self.request.user
                    s.save()
                for s in states.deleted_objects:
                    if self.object.default_state == s:
                        self.object.default_state = None
                    s.delete()

                lineage.instance = self.object
                for lin in lineage.save(commit=False):
                    # if the form has been cleared and there are no children,
                    # let's clean up a little
                    if not (lin.mutation and lin.parent) and not lin.children.exists():
                        lin.delete()
                    else:
                        if not lin.created_by:
                            lin.created_by = self.request.user
                        lin.updated_by = self.request.user
                        lin.reference = self.object.primary_reference
                        lin.save()
                for lin in lineage.deleted_objects:
                    lin.delete()

                if hasattr(self.object, "lineage"):
                    if not self.object.seq:
                        seq = self.object.lineage.parent.protein.seq.mutate(self.object.lineage.mutation)
                        self.object.seq = str(seq)
                    if not self.object.parent_organism:
                        self.object.parent_organism = self.object.lineage.root_node.protein.parent_organism

                comment = f"{self.object} {self.get_form_type()} form."
                chg_string = "\n".join(get_form_changes(form, states, lineage))

                if not self.request.user.is_staff:
                    self.object.status = "pending"
                    msg = f"User: {self.request.user.username}\n"
                    f"Protein: {self.object}\n\n{chg_string}\n\n"
                    f"{self.request.build_absolute_uri(self.object.get_absolute_url())}"
                    mail_managers(comment, msg, fail_silently=True)
                # else:
                #     self.object.status = 'approved'

                self.object.save()
                reversion.set_user(self.request.user)
                reversion.set_comment(chg_string)
                try:
                    uncache_protein_page(self.object.slug, self.request)
                except Exception as e:
                    logger.error("failed to uncache protein: %s", e)

                check_switch_type(self.object, self.request)

        return HttpResponseRedirect(self.get_success_url())


# https://medium.com/@adandan01/django-inline-formsets-example-mybook-420cc4b6225d
class ProteinCreateView(ProteinCreateUpdateMixin, CreateView):
    """renders html for protein submission page"""

    model = Protein
    form_class = ProteinForm

    def get_form(self, *args, **kwargs):
        form = super().get_form(*args, **kwargs)
        # make reference doi required when creating a protein
        form.fields["reference_doi"].required = True
        return form

    def get_context_data(self, **kwargs):
        data = super().get_context_data(**kwargs)
        if self.request.POST:
            data["states"] = StateFormSet(self.request.POST)
            data["lineage"] = LineageFormSet(self.request.POST)
        else:
            data["states"] = StateFormSet()
            data["lineage"] = LineageFormSet()
        return data

    def form_valid(self, form):
        form.instance.created_by = self.request.user
        return super().form_valid(form)


class ProteinUpdateView(ProteinCreateUpdateMixin, UpdateView):
    """renders html for protein submission page"""

    model = Protein
    form_class = ProteinForm

    def get_form(self, form_class=None):
        """Return an instance of the form to be used in this view."""
        if form_class is None:
            form_class = self.get_form_class()
        form_kwargs = self.get_form_kwargs()
        if self.object.primary_reference:
            doi = self.object.primary_reference.doi
            if doi:
                form_kwargs["initial"]["reference_doi"] = doi
        return form_class(**form_kwargs)

    def get_object(self, queryset=None):
        if queryset is None:
            queryset = self.get_queryset()
        try:
            obj = queryset.get(slug=self.kwargs.get("slug"))
        except Protein.DoesNotExist:
            try:
                obj = queryset.get(uuid=self.kwargs.get("slug", "").upper())
            except Protein.DoesNotExist as e:
                raise Http404("No protein found matching this query") from e
        if obj.status == "hidden" and obj.created_by != self.request.user and not self.request.user.is_staff:
            raise Http404("No protein found matching this query")
        return obj

    def get_context_data(self, **kwargs):
        data = super().get_context_data(**kwargs)
        if self.request.POST:
            data["states"] = StateFormSet(self.request.POST, instance=self.object)
            data["states"].full_clean()  # why is this here?

            data["lineage"] = LineageFormSet(self.request.POST, instance=self.object)
        else:
            data["states"] = StateFormSet(instance=self.object)
            data["lineage"] = LineageFormSet(instance=self.object)
        return data

    def form_valid(self, form):
        form.instance.updated_by = self.request.user
        return super().form_valid(form)


class ActivityView(ListView):
    template_name = "proteins/activity.html"
    stateprefetch = Prefetch(
        "states",
        queryset=State.objects.prefetch_related("spectra").order_by("-is_dark", "em_max"),
    )
    queryset = Protein.visible.prefetch_related(stateprefetch, "primary_reference").order_by("-created")[:18]

    def get_context_data(self, **kwargs):
        data = super().get_context_data(**kwargs)
        stateprefetch = Prefetch(
            "states",
            queryset=State.objects.prefetch_related("spectra").order_by("-is_dark", "em_max"),
        )
        data["proteins_by_date"] = (
            Protein.visible.annotate(nstates=Count("states"))
            .prefetch_related(stateprefetch, "primary_reference")
            .order_by(F("primary_reference__date").desc(nulls_last=True))[:15]
        )
        try:
            data["most_viewed"] = {k: v[:12] for k, v in cached_ga_popular().items()}
        except Exception as e:
            logger.error(e)
        data["most_favorited"] = most_favorited(max_results=18)
        return data


@cache_page(60 * 120)
def spectra_image(request, slug, **kwargs):
    protein = get_object_or_404(Protein.objects.select_related("default_state"), slug=slug)
    try:
        d = {}
        for k, v in request.GET.dict().items():
            if k == "xlim":
                tmp = [int(x) for x in v.split(",")][:2]
                if len(tmp) == 1:
                    tmp.append(750)
                d[k] = tmp
            elif v.lower() in ("false", "no"):
                d[k] = False
            elif v.lower() in ("true", "yes"):
                d[k] = True
            else:
                d[k] = int(v) if v.isdigit() else v
    except Exception:
        return HttpResponseBadRequest("failed to parse url parameters as JSON")
    try:
        fmt = kwargs.get("extension", "png")
        byt = protein.spectra_img(fmt, output=io.BytesIO(), **d)
    except Exception as e:
        logger.error(e)
        return HttpResponseBadRequest(f"failed to parse url parameters as JSON: {e}")
    if byt:
        byt.seek(0)
        if fmt == "svg":
            fmt += "+xml"
        return HttpResponse(byt, content_type=f"image/{fmt}")
    raise Http404()


def protein_table(request):
    """Renders html for protein table page.

    The table is now a React app that fetches data from the API endpoint.
    """
    return render(request, "table.html")


class ComparisonView(base.TemplateView):
    template_name = "compare.html"

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        if "proteins" in self.kwargs:
            # the page was requested with slugs in the URL, maintain that order
            ids = kwargs.get("proteins", "").split(",")
            p = Case(*[When(slug=slug, then=pos) for pos, slug in enumerate(ids)])
            prots = Protein.objects.filter(slug__in=ids).prefetch_related("states__spectra").order_by(p)
        else:
            # try to use chronological order
            ids = self.request.session.get("comparison", [])
            prots = (
                Protein.objects.filter(slug__in=ids)
                .prefetch_related("states__spectra")
                .order_by("primary_reference__year")
            )
        prots_with_seqs = prots.exclude(seq__isnull=True)
        if prots_with_seqs.count() > 2:
            a, _ = prots_with_seqs.to_tree("html")
            alignment = "\n".join(a.splitlines()[3:]).replace("FFEEE0", "FFFFFF")
            for p in prots:
                alignment = alignment.replace(f"{p.uuid:16}", f"{p.name[:16]:16}")
            context["alignment"] = alignment
        elif prots_with_seqs.count() == 2:
            seq_a = prots_with_seqs[0]
            seq_b = prots_with_seqs[1]
            if seq_a.seq and seq_b.seq:
                algn = seq_a.seq.align_to(seq_b.seq)
                out = []
                for i, row in enumerate(str(algn).splitlines()):
                    head = ""
                    if i % 4 == 0:
                        head = prots[0].name
                    elif i % 2 == 0:
                        head = prots[1].name
                    out.append("{:<18.16}{}".format(head if len(head) < 17 else head[:13] + "...", row))
                out.append("\n")
                context["alignment"] = "\n".join(out)
                context["mutations"] = str(seq_a.seq.mutations_to(seq_b.seq))
        else:
            context["alignment"] = None
        prots = list(prots)
        context["proteins"] = prots
        refs = Reference.objects.filter(primary_proteins__in=prots)
        refs = refs | Reference.objects.filter(proteins__in=prots)
        refs = refs.distinct("id").order_by("id")
        context["references"] = sorted(refs, key=lambda x: x.year)
        spectra = []
        for prot in prots:
            for state in prot.states.all():
                spectra.extend(state.spectra.all())
        context["spectra_ids"] = ",".join([str(sp.id) for sp in spectra])
        context["hidden_spectra"] = ",".join([str(sp.id) for sp in spectra if sp.subtype in ("2p")])

        return context


def problems_gaps(request):
    return render(
        request,
        "problems_gaps.html",
        {
            "noseqs": Protein.objects.filter(seq__isnull=True).values("name", "slug"),
            "nostates": Protein.objects.filter(states=None).values("name", "slug"),
            "noparent": Protein.objects.filter(parent_organism__isnull=True),
            "only2p": (
                State.objects.filter(spectra__subtype="2p")
                .exclude(spectra__subtype="ex")
                .distinct("protein")
                .values("protein__name", "protein__slug")
            ),
            "nolineage": Protein.objects.filter(lineage=None).annotate(ns=Count("states__spectra")).order_by("-ns"),
            "request": request,
        },
    )


def problems_inconsistencies(request):
    titles = reduce(
        operator.or_,
        (Q(primary_reference__title__icontains=item) for item in ["activat", "switch", "convert", "dark", "revers"]),
    )
    names = reduce(
        operator.or_,
        (Q(name__startswith=item) for item in ["PA", "rs", "mPA", "PS", "mPS"]),
    )
    switchers = Protein.objects.exclude(name__startswith="mCerulean").filter(titles)
    switchers = switchers | Protein.objects.filter(names)
    switchers = switchers.annotate(ns=Count("states")).filter(ns=1)

    gb_mismatch = []
    with_genbank = Protein.objects.exclude(genbank=None).exclude(seq=None).values("slug", "name", "genbank", "seq")
    gbseqs = get_cached_gbseqs([g["genbank"] for g in with_genbank])
    for item in with_genbank:
        if item["genbank"] in gbseqs:
            gbseq = gbseqs[item["genbank"]][0]
            ourseq = item["seq"]
            if ourseq != gbseq:
                gb_mismatch.append(
                    (
                        item["name"],
                        item["slug"],
                        item["genbank"],
                        str(ourseq.mutations_to(gbseq)),
                        ourseq,
                    )
                )
    p = list(
        Protein.objects.annotate(ndark=Count("states", filter=Q(states__is_dark=True)))
        .annotate(nfrom=Count("transitions__from_state", distinct=True))
        .prefetch_related("states", "transitions")
    )
    bad_switch = []
    for prot in p:
        suggestion = suggested_switch_type(prot)
        if (prot.switch_type or suggestion) and (prot.switch_type != suggestion):
            bad_switch.append((prot, dict(Protein.SWITCHING_CHOICES).get(suggestion)))

    return render(
        request,
        "problems_inconsistencies.html",
        {
            "histags": Protein.objects.filter(seq__icontains="HHHHH").values("name", "slug"),
            "linprobs": [(node.protein, v) for node, v in check_lineages()[0].items()],
            "nomet": Protein.objects.exclude(seq__isnull=True).exclude(seq__istartswith="M"),
            "bad_switch": bad_switch,
            "switchers": switchers,
            "request": request,
            "mismatch": gb_mismatch,
        },
    )


def add_reference(request, slug=None):
    try:
        with reversion.create_revision():
            doi = request.POST.get("reference_doi").lower()
            p = Protein.objects.get(slug=slug)
            ref, _created = Reference.objects.get_or_create(doi=doi)
            p.references.add(ref)
            if not request.user.is_staff:
                p.status = "pending"
                msg = (
                    f"User: {request.user.username}\nProtein: {p}\nReference: {ref}. {ref.title}\n\n"
                    f"{request.build_absolute_uri(p.get_absolute_url())}"
                )
                mail_managers("Reference Added", msg, fail_silently=True)
            p.save()
            reversion.set_user(request.user)
            reversion.set_comment(f"Added Reference: {ref}")
            try:
                uncache_protein_page(slug, request)
            except Exception as e:
                logger.error("failed to uncache protein: %s", e)
        return JsonResponse({"status": "success"})
    except Exception as e:
        return JsonResponse({"status": "failed", "msg": e})


def add_protein_excerpt(request, slug=None):
    try:
        with reversion.create_revision():
            doi = request.POST.get("excerpt_doi").lower()
            p = Protein.objects.get(slug=slug)
            content = request.POST.get("excerpt_content")
            if content:
                ref, _created = Reference.objects.get_or_create(doi=doi)
                p.references.add(ref)
                excerpt = Excerpt.objects.create(reference=ref, content=strip_tags(content), created_by=request.user)
                excerpt.proteins.add(p)
                if not request.user.is_staff:
                    msg = (
                        f"User: {request.user.username}\nProtein: {p}\nReference: {ref}, {ref.title}\nExcerpt: "
                        f"{strip_tags(content)}\n{request.build_absolute_uri(p.get_absolute_url())}"
                    )
                    mail_managers("Excerpt Added", msg, fail_silently=True)
                p.save()
                reversion.set_user(request.user)
                reversion.set_comment(f"Added Excerpt from {ref}")
                try:
                    uncache_protein_page(slug, request)
                except Exception as e:
                    logger.error("failed to uncache protein: %s", e)
        return JsonResponse({"status": "success"})
    except Exception as e:
        return JsonResponse({"status": "failed", "msg": e})


@staff_member_required
def revert_version(request, ver=None):
    with contextlib.suppress(Exception):
        version = Version.objects.get(id=ver)
        version.revision.revert(delete=True)
        return JsonResponse({})


@staff_member_required
def revert_revision(request, rev=None):
    revision = get_object_or_404(Revision, id=rev)

    with transaction.atomic():
        revision.revert(delete=True)
        proteins = {v.object for v in revision.version_set.all() if v.object._meta.model_name == "protein"}
        if len(proteins) == 1:
            p = proteins.pop()
            with reversion.create_revision():
                reversion.set_user(request.user)
                reversion.set_comment(f"Reverted to revision dated {revision.date_created}")
                p.save()
                with contextlib.suppress(Exception):
                    uncache_protein_page(p.slug, request)

    return JsonResponse({"status": 200})


@login_required
def update_transitions(request, slug=None):
    template_name = "proteins/forms/_transition_form.html"
    obj = Protein.objects.get(slug=slug)
    if request.method == "POST":
        formset = StateTransitionFormSet(request.POST, instance=obj)
        if not formset.is_valid():
            return render(request, template_name, {"transition_form": formset}, status=422)

        with transaction.atomic():
            with reversion.create_revision():
                formset.save()
                chg_string = "\n".join(get_form_changes(formset))

                if not request.user.is_staff:
                    obj.status = "pending"
                    mail_managers(
                        "Transition updated",
                        f"User: {request.user.username}\nProtein: {obj}\n\n{chg_string}",
                        fail_silently=True,
                    )
                obj.save()
                reversion.set_user(request.user)
                reversion.set_comment(chg_string)
                try:
                    uncache_protein_page(slug, request)
                except Exception as e:
                    logger.error("failed to uncache protein: %s", e)
                check_switch_type(obj, request)
        return HttpResponse(status=200)
    else:
        formset = StateTransitionFormSet(instance=obj)
        formset.form.base_fields["from_state"].queryset = State.objects.filter(protein=obj)
        formset.form.base_fields["to_state"].queryset = State.objects.filter(protein=obj)
        return render(request, template_name, {"transition_form": formset})


@login_required
def protein_bleach_formsets(request, slug):
    template_name = "proteins/protein_bleach_form.html"
    BleachMeasurementFormSet = modelformset_factory(BleachMeasurement, BleachMeasurementForm, extra=1, can_delete=True)
    protein = get_object_or_404(Protein, slug=slug)
    qs = BleachMeasurement.objects.filter(state__protein=protein)
    if request.method == "POST":
        formset = BleachMeasurementFormSet(request.POST, queryset=qs)
        formset.form.base_fields["state"].queryset = State.objects.filter(protein__slug=slug)
        if not formset.is_valid():
            return render(request, template_name, {"formset": formset, "protein": protein})

        with transaction.atomic():
            with reversion.create_revision():
                saved = formset.save(commit=False)
                for s in saved:
                    if not s.created_by:
                        s.created_by = request.user
                    s.updated_by = request.user
                    s.save()
                for s in formset.deleted_objects:
                    s.delete()

                chg_string = "\n".join(get_form_changes(formset))

                if not request.user.is_staff:
                    protein.status = "pending"
                    mail_managers(
                        "BleachMeasurement Added",
                        f"User: {request.user.username}\nProtein: {protein}\n{chg_string}\n\n"
                        f"{request.build_absolute_uri(protein.get_absolute_url())}",
                        fail_silently=True,
                    )
                # else:
                #     protein.status = 'approved'

                protein.save()
                reversion.set_user(request.user)
                reversion.set_comment(chg_string)
                try:
                    uncache_protein_page(slug, request)
                except Exception as e:
                    logger.error("failed to uncache protein: %s", e)
        return HttpResponseRedirect(protein.get_absolute_url())
    else:
        formset = BleachMeasurementFormSet(queryset=qs)
        formset.form.base_fields["state"].queryset = State.objects.filter(protein__slug=slug)
    return render(request, template_name, {"formset": formset, "protein": protein})


@login_required
def bleach_comparison(request, pk=None):
    template_name = "proteins/bleach_comparison_form.html"
    if request.method == "POST":
        formset = bleach_items_formset(request.POST)
        bcf = BleachComparisonForm(request.POST)
        if not formset.is_valid() or not bcf.is_valid():
            return render(request, template_name, {"formset": formset, "mainform": bcf})
        d = bcf.cleaned_data
        d["reference"], _ = Reference.objects.get_or_create(doi=d.pop("reference_doi"))
        for form in formset.forms:
            if form.has_changed():
                d["state"] = form.cleaned_data["state"]
                d["rate"] = form.cleaned_data["rate"]
                BleachMeasurement.objects.create(**d)
        return redirect(d["reference"])
    else:
        formset = bleach_items_formset()
        if pk:
            reference = get_object_or_404(Reference, id=pk)
            bcf = BleachComparisonForm({"reference_doi": reference.doi})
            bcf.fields["reference_doi"].widget.attrs["readonly"] = True
        else:
            bcf = BleachComparisonForm()
    return render(request, template_name, {"formset": formset, "mainform": bcf})


class OrganismListView(ListView):
    """renders html for single reference page"""

    queryset = Organism.objects.annotate(num_prot=Count("proteins"))


class OrganismDetailView(DetailView):
    """renders html for single reference page"""

    queryset = Organism.objects.all().prefetch_related("proteins__states")


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


@login_required
def flag_object(request):
    if not is_ajax(request):
        return HttpResponseNotAllowed([])
    try:
        model_type = request.POST["target_model"]
        id_ = request.POST["target_id"]
        model = apps.get_model(model_type)
        obj = get_object_or_404(model, id=id_)

        status = None
        if request.POST["flagged"] == "1":
            if request.user.is_staff:
                obj.status = model.STATUS.approved
                status = "unflagged"
        else:
            obj.status = model.STATUS.flagged
            status = "flagged"
        obj.save()

        if status:
            try:
                uncache_protein_page(obj.protein.slug, request)
            except Exception as e:
                logger.error("failed to uncache protein: %s", e)

            mail_admins(
                f"{model_type} {status}",
                f"User: {request.user.username}\nObject: {obj}\nID: {obj.pk}\n"
                f"{request.build_absolute_uri(obj.get_absolute_url())}",
                fail_silently=True,
            )

        return JsonResponse({"status": status})
    except (KeyError, ValueError):
        return HttpResponseBadRequest()


def protein_history(request, slug):
    protein = get_object_or_404(Protein, slug=slug)
    return render(
        request,
        "history.html",
        {
            "protein": protein,
            "history": protein.history(
                [
                    "modified",
                    "created_by_id",
                    "status_changed",
                    "emhex",
                    "exhex",
                    "updated_by_id",
                    "status",
                    "seq_comment",
                ]
            ),
            "request": request,
        },
    )


def _get_protein_changes(protein):
    """
    Compare a pending protein with its last approved version and return structured changes.

    Returns a dict with categories of changes:
    - protein_fields: Direct protein field changes
    - states: State additions/removals/modifications
    - transitions: StateTransition changes
    - oser: OSER measurement changes
    - lineage: Lineage changes
    - references: Reference changes
    """
    from collections import defaultdict
    from reversion.models import Version

    changes = {
        "protein_fields": {},
        "states": {"added": [], "removed": [], "modified": {}},
        "spectra": {"added": [], "removed": [], "modified": {}},
        "transitions": {"added": [], "removed": [], "modified": {}},
        "oser": {"added": [], "removed": [], "modified": {}},
        "lineage": {},
        "references": {"added": [], "removed": []},
        "bleach": {"added": [], "removed": [], "modified": {}},
    }

    # Get last approved version
    last_approved_version = protein.last_approved_version()
    if not last_approved_version:
        # No approved version - this is a new protein submission
        changes["is_new"] = True
        return changes

    changes["is_new"] = False

    # If last_approved is a Version object, we need to restore the protein from it
    if isinstance(last_approved_version, Version):
        # Use the old_object helper from history.py to restore the protein
        from proteins.util.history import old_object

        last_approved = old_object(last_approved_version)
    else:
        # It's already a Protein instance (when status == "approved")
        last_approved = last_approved_version

    # Compare protein-level fields
    field_map = {
        "name": "Name",
        "seq": "Sequence",
        "chromophore": "Chromophore",
        "aliases": "Aliases",
        "pdb": "PDB IDs",
        "genbank": "GenBank",
        "uniprot": "UniProt",
        "ipg_id": "IPG ID",
        "mw": "Molecular Weight",
        "agg": "Oligomerization",
        "oser": "OSER",
        "switch_type": "Switching Type",
        "blurb": "Description",
        "cofactor": "Cofactor",
        "seq_validated": "Sequence Validated",
        "seq_comment": "Sequence Comment",
        "parent_organism_id": "Parent Organism",
        "primary_reference_id": "Primary Reference",
        "default_state_id": "Default State",
    }

    for field, label in field_map.items():
        old_val = getattr(last_approved, field, None)
        new_val = getattr(protein, field, None)

        # Special handling for different field types
        if field == "seq":
            if old_val != new_val:
                if old_val and new_val:
                    # Show mutations
                    mutations = old_val.mutations_to(new_val) if hasattr(old_val, "mutations_to") else None
                    changes["protein_fields"][label] = {
                        "old": str(old_val)[:50] + "..." if len(str(old_val)) > 50 else str(old_val),
                        "new": str(new_val)[:50] + "..." if len(str(new_val)) > 50 else str(new_val),
                        "mutations": str(mutations) if mutations else None,
                    }
                else:
                    changes["protein_fields"][label] = {
                        "old": str(old_val) if old_val else None,
                        "new": str(new_val) if new_val else None,
                    }
        elif field == "agg":
            if old_val != new_val:
                old_display = last_approved.get_agg_display() if old_val else None
                new_display = protein.get_agg_display() if new_val else None
                if old_display != new_display:
                    changes["protein_fields"][label] = {"old": old_display, "new": new_display}
        elif field == "switch_type":
            if old_val != new_val:
                old_display = last_approved.get_switch_type_display() if old_val else None
                new_display = protein.get_switch_type_display() if new_val else None
                if old_display != new_display:
                    changes["protein_fields"][label] = {"old": old_display, "new": new_display}
        elif field == "cofactor":
            if old_val != new_val:
                old_display = last_approved.get_cofactor_display() if old_val else None
                new_display = protein.get_cofactor_display() if new_val else None
                if old_display != new_display:
                    changes["protein_fields"][label] = {"old": old_display, "new": new_display}
        elif field == "parent_organism_id":
            if old_val != new_val:
                old_org = last_approved.parent_organism if old_val else None
                new_org = protein.parent_organism if new_val else None
                changes["protein_fields"][label] = {
                    "old": str(old_org) if old_org else None,
                    "new": str(new_org) if new_org else None,
                }
        elif field in ["primary_reference_id", "default_state_id"]:
            # Skip these for now - we'll handle them separately if needed
            continue
        else:
            if old_val != new_val:
                changes["protein_fields"][label] = {"old": old_val, "new": new_val}

    # Compare states
    old_states = {s.id: s for s in last_approved.states.all()}
    new_states = {s.id: s for s in protein.states.all()}

    # Find added states (states that exist now but didn't before)
    for state_id, state in new_states.items():
        if state_id not in old_states:
            # Check if this is truly new or just has a new ID
            # Match by unique fields: name, ex_max, em_max, ext_coeff, qy
            is_truly_new = True
            for old_state in old_states.values():
                if (
                    old_state.name == state.name
                    and old_state.ex_max == state.ex_max
                    and old_state.em_max == state.em_max
                    and old_state.ext_coeff == state.ext_coeff
                    and old_state.qy == state.qy
                ):
                    is_truly_new = False
                    break

            if is_truly_new:
                changes["states"]["added"].append(
                    {
                        "name": str(state),
                        "ex_max": state.ex_max,
                        "em_max": state.em_max,
                        "ext_coeff": state.ext_coeff,
                        "qy": state.qy,
                    }
                )

    # Find removed states
    for state_id, old_state in old_states.items():
        if state_id not in new_states:
            # Check if truly removed or just has new ID
            is_truly_removed = True
            for state in new_states.values():
                if (
                    old_state.name == state.name
                    and old_state.ex_max == state.ex_max
                    and old_state.em_max == state.em_max
                    and old_state.ext_coeff == state.ext_coeff
                    and old_state.qy == state.qy
                ):
                    is_truly_removed = False
                    break

            if is_truly_removed:
                changes["states"]["removed"].append(
                    {
                        "name": str(old_state),
                        "ex_max": old_state.ex_max,
                        "em_max": old_state.em_max,
                    }
                )

    # Compare modified states (match by ID or by unique fields)
    state_field_map = {
        "name": "Name",
        "ex_max": "Ex max",
        "em_max": "Em max",
        "ext_coeff": "Extinction Coefficient",
        "qy": "Quantum Yield",
        "pka": "pKa",
        "maturation": "Maturation",
        "lifetime": "Lifetime",
        "is_dark": "Is Dark",
        "twop_ex_max": "2P Ex max",
        "twop_peakGM": "2P Peak GM",
        "twop_qy": "2P QY",
    }

    for state_id in set(old_states.keys()) & set(new_states.keys()):
        old_state = old_states[state_id]
        new_state = new_states[state_id]
        state_changes = {}

        for field, label in state_field_map.items():
            old_val = getattr(old_state, field, None)
            new_val = getattr(new_state, field, None)
            if old_val != new_val:
                state_changes[label] = {"old": old_val, "new": new_val}

        if state_changes:
            changes["states"]["modified"][str(new_state)] = state_changes

    # Compare references
    old_refs = set(last_approved.references.values_list("id", flat=True))
    new_refs = set(protein.references.values_list("id", flat=True))

    added_refs = new_refs - old_refs
    removed_refs = old_refs - new_refs

    if added_refs:
        from references.models import Reference

        for ref_id in added_refs:
            ref = Reference.objects.get(id=ref_id)
            changes["references"]["added"].append(str(ref.citation))

    if removed_refs:
        from references.models import Reference

        for ref_id in removed_refs:
            ref = Reference.objects.get(id=ref_id)
            changes["references"]["removed"].append(str(ref.citation))

    # Compare lineage
    try:
        old_lineage = last_approved.lineage
        new_lineage = protein.lineage
        lineage_fields = {"parent": "Parent", "mutation": "Mutation", "reference": "Reference"}

        for field, label in lineage_fields.items():
            old_val = getattr(old_lineage, field, None)
            new_val = getattr(new_lineage, field, None)
            if old_val != new_val:
                changes["lineage"][label] = {
                    "old": str(old_val) if old_val else None,
                    "new": str(new_val) if new_val else None,
                }
    except Exception:
        # No lineage or error comparing
        pass

    # Compare OSER measurements
    old_osers = {o.id: o for o in last_approved.oser_measurements.all()}
    new_osers = {o.id: o for o in protein.oser_measurements.all()}

    for oser_id in set(new_osers.keys()) - set(old_osers.keys()):
        oser = new_osers[oser_id]
        changes["oser"]["added"].append(str(oser))

    for oser_id in set(old_osers.keys()) - set(new_osers.keys()):
        oser = old_osers[oser_id]
        changes["oser"]["removed"].append(str(oser))

    # Compare transitions
    old_transitions = {t.id: t for t in last_approved.transitions.all()}
    new_transitions = {t.id: t for t in protein.transitions.all()}

    for trans_id in set(new_transitions.keys()) - set(old_transitions.keys()):
        trans = new_transitions[trans_id]
        changes["transitions"]["added"].append(str(trans))

    for trans_id in set(old_transitions.keys()) - set(new_transitions.keys()):
        trans = old_transitions[trans_id]
        changes["transitions"]["removed"].append(str(trans))

    # Clean up empty sections
    changes = {k: v for k, v in changes.items() if v and (not isinstance(v, dict) or any(v.values()))}

    return changes


@staff_member_required
def pending_proteins_dashboard(request):
    """Dashboard for reviewing pending protein submissions and changes."""
    from django.contrib.auth.decorators import permission_required

    # Get all pending proteins
    pending_proteins = (
        Protein.objects.filter(status="pending")
        .select_related(
            "created_by",
            "updated_by",
            "parent_organism",
            "primary_reference",
            "default_state",
        )
        .prefetch_related(
            "states",
            "states__spectra",
            "states__bleach_measurements",
            "transitions",
            "oser_measurements",
            "references",
        )
        .order_by("-modified")
    )

    proteins_data = []
    for protein in pending_proteins:
        # Get changes for this protein
        changes = _get_protein_changes(protein)

        # Only include proteins that have changes (or are new)
        if changes.get("is_new") or any(v for k, v in changes.items() if k != "is_new"):
            proteins_data.append(
                {
                    "id": protein.id,
                    "slug": protein.slug,
                    "name": protein.name,
                    "uuid": protein.uuid,
                    "created": protein.created,
                    "modified": protein.modified,
                    "created_by": protein.created_by,
                    "updated_by": protein.updated_by,
                    "created_by_email": protein.created_by.email if protein.created_by else None,
                    "is_new": changes.get("is_new", False),
                    "changes": changes,
                    "admin_url": f"/admin/proteins/protein/{protein.id}/change/",
                    "detail_url": protein.get_absolute_url(),
                }
            )

    context = {
        "proteins": proteins_data,
        "count": len(proteins_data),
    }

    return render(request, "pending_proteins_dashboard.html", context)


@staff_member_required
@reversion.create_revision()
def pending_protein_action(request):
    """Handle actions (approve/reject) on pending proteins."""
    from django.views.decorators.http import require_POST
    from django.http import JsonResponse

    if request.method != "POST":
        return JsonResponse({"success": False, "error": "POST required"}, status=405)

    try:
        protein_ids = request.POST.getlist("protein_ids[]")
        action = request.POST.get("action")

        if not protein_ids or not action:
            return JsonResponse({"success": False, "error": "Missing protein_ids or action"}, status=400)

        proteins = Protein.objects.filter(id__in=protein_ids, status="pending")

        if not proteins.exists():
            return JsonResponse({"success": False, "error": "No pending proteins found with provided IDs"}, status=404)

        count = proteins.count()

        if action == "approve":
            # Set revision comment
            reversion.set_comment(f"Approved {count} protein(s) via moderation dashboard")
            reversion.set_user(request.user)

            proteins.update(status="approved")
            # Clear cache for affected protein pages
            for protein in proteins:
                with contextlib.suppress(Exception):
                    uncache_protein_page(protein.slug, request)
            message = f"Approved {count} protein(s)"

        elif action == "reject":
            # Set revision comment
            reversion.set_comment(f"Rejected {count} protein(s) via moderation dashboard")
            reversion.set_user(request.user)

            # For reject, we might want to revert to last approved version
            # or just mark as hidden. Let's mark as hidden for now.
            proteins.update(status="hidden")
            message = f"Rejected (hidden) {count} protein(s)"

        else:
            return JsonResponse({"success": False, "error": f"Unknown action: {action}"}, status=400)

        return JsonResponse({"success": True, "message": message, "count": count})

    except Exception as e:
        logger.exception("Error in pending_protein_action: %s", e)
        return JsonResponse({"success": False, "error": str(e)}, status=500)
