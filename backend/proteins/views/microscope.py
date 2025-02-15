import json
from collections import defaultdict
from urllib.parse import quote

from django.contrib import messages
from django.contrib.auth import get_user_model
from django.contrib.contenttypes.models import ContentType
from django.contrib.messages.views import SuccessMessageMixin
from django.contrib.postgres.aggregates import ArrayAgg
from django.core.exceptions import PermissionDenied
from django.core.mail import mail_admins
from django.db import transaction
from django.db.models import CharField, Count, F, Q, Value
from django.db.models.functions import Lower
from django.http import (
    Http404,
    HttpResponseNotAllowed,
    HttpResponseRedirect,
    JsonResponse,
)
from django.urls import resolve, reverse_lazy
from django.utils.decorators import method_decorator
from django.views.decorators.clickjacking import xframe_options_exempt
from django.views.generic import (
    CreateView,
    DeleteView,
    DetailView,
    ListView,
    UpdateView,
)

from fpbase.celery import app
from fpbase.util import is_ajax

from ..forms import MicroscopeForm, OpticalConfigFormSet
from ..models import (
    Camera,
    Dye,
    Light,
    Microscope,
    OpticalConfig,
    ProteinCollection,
    Spectrum,
    State,
)

# from ..util.efficiency import microscope_efficiency_report
from ..models.efficiency import OcFluorEff
from ..tasks import calculate_scope_report
from .mixins import OwnableObject


def update_scope_report(request):
    job_id = request.POST.get("job_id")
    scope_id = request.POST.get("scope_id")

    if request.POST.get("action") == "update":
        outdated = request.POST.get("outdated")
        if outdated:
            outdated = json.loads(outdated)
        if scope_id:
            try:
                # this is throwing connection resets
                active = app.control.inspect().active()
            except Exception:
                active = None
            if active:
                for _worker, jobs in active.items():
                    for job in jobs:
                        if job["name"].endswith("calculate_scope_report") and (scope_id in job["args"]):
                            return JsonResponse({"status": 200, "job": job["id"]})
                    if len(jobs) >= 4:
                        return JsonResponse({"status": 200, "job": None, "waiting": True})
            job_id = calculate_scope_report.delay(scope_id, outdated_ids=outdated).id
            return JsonResponse({"status": 200, "job": job_id})
    elif request.POST.get("action") == "check":
        if job_id:
            result = app.AsyncResult(job_id)
            ready = result.ready()
            return JsonResponse({"status": 200, "ready": ready, "info": result.info})
    elif request.POST.get("action") == "cancel":
        if job_id:
            result = app.AsyncResult(job_id)
            result.revoke(terminate=True)
            return JsonResponse(
                {
                    "status": 200,
                    "ready": result.ready(),
                    "info": "result.info",
                    "canceled": True,
                }
            )
    return JsonResponse({"status": 404})


def scope_report_json(request, pk):
    microscope = Microscope.objects.get(id=pk)
    oclist = microscope.optical_configs.values_list("id")

    state_ct = ContentType.objects.get_for_model(State)
    dye_ct = ContentType.objects.get_for_model(Dye)

    effs = list(
        OcFluorEff.objects.exclude(ex_eff=None)
        .filter(content_type=state_ct, oc__in=oclist)
        .annotate(
            fluor_id=F("state__protein__uuid"),
            fluor_slug=F("state__slug"),
            type=Value("p", CharField()),
        )
        .values(
            "fluor_id",
            "fluor_name",
            "fluor_slug",
            "ex_eff",
            "em_eff",
            "ex_eff_broad",
            "brightness",
            "type",
            "oc__name",
        )
    )

    effs.extend(
        list(
            OcFluorEff.objects.exclude(ex_eff=None)
            .filter(content_type=dye_ct, oc__in=oclist)
            .annotate(
                fluor_id=F("dye__id"),
                fluor_slug=F("dye__slug"),
                type=Value("d", CharField()),
            )
            .values(
                "fluor_id",
                "fluor_name",
                "fluor_slug",
                "ex_eff",
                "em_eff",
                "ex_eff_broad",
                "brightness",
                "type",
                "oc__name",
            )
        )
    )

    data = defaultdict(list)
    for item in effs:
        if not (item["ex_eff"] and item["em_eff"]):
            continue
        data[item["oc__name"]].append(
            {
                "fluor": item["fluor_name"],
                "fluor_slug": item["fluor_slug"],
                "fluor_id": item["fluor_id"],
                "ex_eff": item["ex_eff"],
                "ex_eff_broad": item["ex_eff_broad"],
                "em_eff": item["em_eff"],
                "brightness": item["brightness"] or None,
                "shape": "circle" if item["type"] == "p" else "square",
                "url": microscope.get_absolute_url()
                + "?c={}&p={}".format(quote(item["oc__name"]), quote(item["fluor_slug"])),
            }
        )

    states = State.objects.with_spectra().prefetch_related("protein")
    fluors = {
        i["slug"]: i
        for i in states.annotate(
            agg=F("protein__agg"),
            type=Value("p", CharField()),
            switch_type=F("protein__switch_type"),
            uuid=F("protein__uuid"),
        ).values(
            "slug",
            "ext_coeff",
            "qy",
            "agg",
            "emhex",
            "ex_max",
            "em_max",
            "type",
            "switch_type",
            "uuid",
        )
    }
    dyes = Dye.objects.with_spectra()
    fluors.update(
        {
            i["slug"]: i
            for i in dyes.annotate(
                uuid=F("id"),
                agg=Value("", CharField()),
                type=Value("d", CharField()),
                switch_type=Value("", CharField()),
            ).values(
                "slug",
                "ext_coeff",
                "qy",
                "agg",
                "emhex",
                "ex_max",
                "em_max",
                "type",
                "switch_type",
                "uuid",
            )
        }
    )
    return JsonResponse(
        {
            "microscope": pk,
            "report": [{"key": k, "values": v} for k, v in data.items()],
            "fluors": fluors,
        }
    )


class ScopeReportView(DetailView):
    template_name = "proteins/scope_report.html"
    queryset = Microscope.objects.all()

    def post(self, request, *args, **kwargs):
        if is_ajax(request):
            return update_scope_report(request)
        return HttpResponseNotAllowed([])

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        probe_count = State.objects.with_spectra().count() + Dye.objects.with_spectra().count()
        ids = self.object.optical_configs.all().values_list("id", flat=True)
        effs = OcFluorEff.objects.filter(oc__in=ids)
        context["outdated"] = list(effs.outdated().values_list("id", flat=True))
        context["needs_update"] = bool(context["outdated"]) or (probe_count * len(ids) > effs.count())
        cols = (
            ProteinCollection.objects.filter(Q(private=False) | Q(owner_id=self.request.user.id))
            .exclude(proteins=None)
            .annotate(uuids=ArrayAgg("proteins__uuid"))
            .order_by(Lower("name"))
            .values("id", "name", "uuids")
        )
        for c in cols:
            c["uuids"] = json.dumps(c["uuids"])
        context["collections"] = cols
        return context


class MicroscopeCreateUpdateMixin:
    def get_form_type(self):
        return self.request.resolver_match.url_name

    def form_valid(self, form):
        context = self.get_context_data()
        ocformset = context["optical_configs"]

        if not ocformset.is_valid():
            messages.add_message(
                self.request,
                messages.ERROR,
                "Please correct errors in the individual configs tab",
            )
            context.update({"optical_configs": ocformset, "formsets_invalid": "true"})
            return self.render_to_response(context)

        # enforce at least one valid optical config
        ocform_has_forms = any(f.cleaned_data.get("name") for f in ocformset.forms if not f.cleaned_data.get("DELETE"))
        if not (ocform_has_forms or form.cleaned_data.get("optical_configs")):
            messages.add_message(
                self.request,
                messages.ERROR,
                "What good is a microscope without any optical "
                "configurations? Please add at least one config "
                "on either tab below",
            )
            context.update({"optical_configs": ocformset})
            return self.render_to_response(context)

        self.object = form.save(commit=False)

        # otherwise save...
        with transaction.atomic():
            ocformset.instance = self.object
            ocformset.save()
            self.object.save()
        return HttpResponseRedirect(self.get_success_url())


class MicroscopeCreateView(MicroscopeCreateUpdateMixin, OwnableObject, CreateView):
    """renders html for microscope creation page"""

    model = Microscope
    form_class = MicroscopeForm

    def form_valid(self, form):
        form.instance.created_by = self.request.user
        response = super().form_valid(form)
        if self.object:
            messages.add_message(
                self.request,
                messages.SUCCESS,
                "Welcome to your new microscope page!  You can find this "
                "page any time in your profile, and share this link with "
                "others.  To <strong>update</strong> or <strong>delete</strong> this microscope, "
                "click the <i class='fas fa-cog mx-1'></i> icon below the graph.",
            )
            if not self.request.user.is_staff:
                try:
                    mail_admins(
                        "Microscope Created",
                        f"User: {self.request.user.username}\nMicroscope: {self.object}"
                        f"\n{self.request.build_absolute_uri(self.object.get_absolute_url())}",
                        fail_silently=True,
                    )
                except Exception:
                    pass
        return response

    def get_context_data(self, **kwargs):
        data = super().get_context_data(**kwargs)
        if self.request.POST:
            data["optical_configs"] = OpticalConfigFormSet(self.request.POST)
        else:
            data["optical_configs"] = OpticalConfigFormSet()
        return data


class MicroscopeUpdateView(SuccessMessageMixin, MicroscopeCreateUpdateMixin, OwnableObject, UpdateView):
    model = Microscope
    form_class = MicroscopeForm
    success_message = "Update successful!"

    def dispatch(self, request, *args, **kwargs):
        if not self.get_object().has_change_permission(self.request):
            raise PermissionDenied
        return super().dispatch(request, *args, **kwargs)

    def form_valid(self, form):
        form.instance.updated_by = self.request.user
        return super().form_valid(form)

    def get_context_data(self, **kwargs):
        data = super().get_context_data(**kwargs)
        if self.request.POST:
            data["optical_configs"] = OpticalConfigFormSet(self.request.POST, instance=self.object)
        else:
            data["optical_configs"] = OpticalConfigFormSet(instance=self.object, queryset=OpticalConfig.objects.all())
        return data


class MicroscopeDetailView(DetailView):
    """renders html for microscope detail/spectrum page"""

    queryset = Microscope.objects.all().prefetch_related(
        "optical_configs__filterplacement_set__filter",
        "optical_configs__camera",
        "optical_configs__light",
    )

    def get_object(self, queryset=None):
        try:
            scope = (
                Microscope.objects.filter(id__istartswith=self.kwargs.get("pk"))
                .prefetch_related(
                    "optical_configs__filterplacement_set__filter",
                    "optical_configs__camera",
                    "optical_configs__light",
                    "optical_configs",
                )
                .get()
            )
            return scope
        except Microscope.MultipleObjectsReturned as e:
            raise Http404("Multiple microscopes found matching this query") from e
        except Microscope.DoesNotExist as e:
            raise Http404("No microscope found matching this query") from e

    def get_context_data(self, **kwargs):
        data = super().get_context_data(**kwargs)
        if not (self.object.cameras.exists() or self.object.extra_cameras.exists()):
            data["cameras"] = Camera.objects.all()
        if not (
            self.object.lights.exists()
            or self.object.lasers
            or self.object.extra_lights.exists()
            or self.object.extra_lasers
        ):
            data["lights"] = Light.objects.all()
        if self.object.collection:
            proteins = self.object.collection.proteins.with_spectra().prefetch_related("states")
            data["probeslugs"] = [
                {"slug": s.slug, "name": str(s)} for p in proteins for s in p.states.all() if s.spectra
            ]
        else:
            data["probeslugs"] = Spectrum.objects.fluorlist()
        if len(self.object.em_filters) + len(self.object.ex_filters) < 20:
            data["scopespectra"] = json.dumps(self.object.spectra_d3())
        else:
            data["scopespectra"] = {}
        # safari = ('Safari' in self.request.META['HTTP_USER_AGENT']
        #           and 'Chrome' not in self.request.META['HTTP_USER_AGENT'])
        # if safari:
        #     messages.add_message(
        #         self.request, messages.INFO, 'Due to slow performance on Safari, '
        #         'wavelength precision has been reduced. '
        #         'Click gear tab for settings, or zoom in for increased performance')
        return data


class MicroscopeEmbedView(MicroscopeDetailView):
    template_name = "proteins/microscope_embed.html"

    @method_decorator(xframe_options_exempt)
    def dispatch(self, *args, **kwargs):
        return super().dispatch(*args, **kwargs)


class MicroscopeList(ListView):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.example_ids = [
            "i6WL2WdgcDMgJYtPrpZcaJ",
            "wKqWbgApvguSNDSRZNSfpN",
            "4yL4ggAozzcMwTU4Ae7zxF",
        ]

    def get_queryset(self):
        # get all collections for current user and all other non-private collections
        qs = Microscope.objects.annotate(nocs=Count("optical_configs")).filter(nocs__gt=1)
        if self.request.user.is_authenticated:
            qs = qs | Microscope.objects.filter(owner=self.request.user)
        if "owner" in self.kwargs:
            qs = qs.filter(owner__username=self.kwargs["owner"])
        else:
            qs = qs.exclude(id__in=self.example_ids)
        return qs.order_by("-created")[:20]

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        if "owner" in self.kwargs:
            owner = self.kwargs["owner"]
            if not get_user_model().objects.filter(username=owner).exists():
                raise Http404()
            context["owner"] = owner
        context["example_list"] = Microscope.objects.filter(id__in=self.example_ids)
        if self.request.user.is_authenticated:
            context["managing"] = Microscope.objects.filter(managers__contains=[self.request.user.email])
        return context


class MicroscopeDeleteView(DeleteView):
    model = Microscope
    success_url = reverse_lazy("proteins:newmicroscope")

    def dispatch(self, request, *args, **kwargs):
        if not self.get_object().owner == self.request.user:
            raise PermissionDenied
        return super().dispatch(request, *args, **kwargs)

    def get_success_url(self):
        redirect_url = reverse_lazy("proteins:microscopes", kwargs={"owner": self.request.user})
        try:
            # check that this is an internal redirection
            resolve(redirect_url)
        except Exception:
            redirect_url = None
        return redirect_url or super().get_success_url()
