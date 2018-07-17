import json
from django.views.generic import TemplateView, DetailView, CreateView, DeleteView, ListView, UpdateView
from django.http import HttpResponseRedirect, Http404
from django.core.mail import mail_admins
from django.urls import reverse_lazy, resolve
from django.db import transaction
from django.contrib.auth import get_user_model
from django.contrib.messages.views import SuccessMessageMixin
from django.utils.decorators import method_decorator
from django.views.decorators.clickjacking import xframe_options_exempt
from django.contrib import messages

from ..models import Microscope, Camera, Light, Spectrum, OpticalConfig, State
from ..forms import MicroscopeForm, OpticalConfigFormSet
from .mixins import OwnableObject
from ..util.efficiency import oclist_efficiency_report


class ScopeReportView(TemplateView):
    template_name = 'proteins/scope_report.html'

    def get_context_data(self, **kwargs):
        context = super(ScopeReportView, self).get_context_data(**kwargs)
        m = Microscope.objects.last()
        from django.core.cache import cache
        if not cache.get('my_key'):
            oc = list(m.optical_configs.all())
            fluors = list(State.objects.with_spectra())
            # fluors.extend(list(Dye.objects.all()))
            report = oclist_efficiency_report(oc, fluors)
            cache.set('my_key', report, 1800)
        _rep = cache.get('my_key')
        report = {}
        data = []
        for oc, dct in _rep.items():
            data.append({
                'key': oc.name,
                'values': []
            })
            for probe, effdata in dct.items():
                if effdata['ex'] and effdata['em']:
                    data[-1]['values'].append({
                        'fluor': probe.fluor_name,
                        'x': effdata['ex'],
                        'y': effdata['em'],
                        'bright': effdata['bright'],
                        'size': effdata['bright'] if effdata['bright'] else 0.01,
                        'color': probe.emhex,
                        'shape': 'circle' if hasattr(probe, 'protein') else 'square',
                        'url': probe.get_absolute_url(),
                    })

        context['report'] = data
        return context


class MicroscopeCreateUpdateMixin:

    def get_form_type(self):
        return self.request.resolver_match.url_name

    def form_valid(self, form):
        context = self.get_context_data()
        ocformset = context['optical_configs']
        self.object = form.save(commit=False)

        if not ocformset.is_valid():
            self.success_message = 'Please correct errors in the individual configs tab'
            context.update({'optical_configs': ocformset, 'formsets_invalid': 'true'})
            return self.render_to_response(context)

        # enforce at least one valid optical config
        ocform_has_forms = any([f.cleaned_data.get('name')
                                for f in ocformset.forms
                                if not f.cleaned_data.get("DELETE")])
        if not (ocform_has_forms or form.cleaned_data.get('optical_configs')):
            self.success_message = ('What good is a microscope without any optical '
                                    'configurations? Please add at least one config '
                                    'on either tab below')
            context.update({'optical_configs': ocformset})
            return self.render_to_response(context)

        # otherwise save...
        with transaction.atomic():
            ocformset.instance = self.object
            ocformset.save()
            if not self.request.user.is_staff:
                mail_admins('Microscope Created',
                            "User: {}\nCollection: {}\n{}".format(
                                self.request.user.username,
                                self.object,
                                self.request.build_absolute_uri(
                                    self.object.get_absolute_url())),
                            fail_silently=True)
            self.object.save()
        return HttpResponseRedirect(self.get_success_url())


class MicroscopeCreateView(SuccessMessageMixin, MicroscopeCreateUpdateMixin,
                           OwnableObject, CreateView):
    ''' renders html for microscope creation page '''
    model = Microscope
    form_class = MicroscopeForm
    success_message = ("Welcome to your new microcope page!  You can find this "
                       "page any time in your profile, and share this link with "
                       "others.  Click the gear icon below the graph to edit "
                       "configurations or delete this microscope.")

    def form_valid(self, form):
        form.instance.created_by = self.request.user
        return super().form_valid(form)

    def get_context_data(self, **kwargs):
        data = super().get_context_data(**kwargs)
        if self.request.POST:
            data['optical_configs'] = OpticalConfigFormSet(self.request.POST)
        else:
            data['optical_configs'] = OpticalConfigFormSet()
        return data


class MicroscopeUpdateView(SuccessMessageMixin, MicroscopeCreateUpdateMixin,
                           OwnableObject, UpdateView):
    model = Microscope
    form_class = MicroscopeForm
    success_message = "Update successful!"

    def form_valid(self, form):
        form.instance.updated_by = self.request.user
        return super().form_valid(form)

    def get_context_data(self, **kwargs):
        data = super().get_context_data(**kwargs)
        if self.request.POST:
            data['optical_configs'] = OpticalConfigFormSet(self.request.POST,
                                                           instance=self.object)
        else:
            data['optical_configs'] = OpticalConfigFormSet(
                instance=self.object,
                queryset=OpticalConfig.objects.all())
        return data


class MicroscopeDetailView(DetailView):
    ''' renders html for microscope detail/spectrum page '''
    queryset = Microscope.objects.all().prefetch_related(
        'optical_configs__filterplacement_set__filter',
        'optical_configs__filters')

    def get_context_data(self, **kwargs):
        data = super().get_context_data(**kwargs)
        if not self.object.cameras.count():
            data['cameras'] = Camera.objects.all()
        if not self.object.lights.count():
            data['lights'] = Light.objects.all()
        if self.object.collection:
            proteins = self.object.collection.proteins.with_spectra() \
                .prefetch_related('proteins', 'proteins__states')
            data['probeslugs'] = [{'slug': s.slug, 'name': str(s)}
                                  for p in proteins for s in p.states.all()
                                  if s.spectra]
        else:
            data['probeslugs'] = Spectrum.objects.fluorlist()
        data['scopespectra'] = json.dumps(self.object.spectra_d3())
        safari = ('Safari' in self.request.META['HTTP_USER_AGENT']
                  and 'Chrome' not in self.request.META['HTTP_USER_AGENT'])
        if safari:
            messages.add_message(
                self.request, messages.INFO, 'Due to slow performance on Safari, '
                'wavelength precision has been reduced. '
                'Click gear tab for settings, or zoom in for increased performance')
        return data


class MicroscopeEmbedView(MicroscopeDetailView):
    template_name = "proteins/microscope_embed.html"

    @method_decorator(xframe_options_exempt)
    def dispatch(self, *args, **kwargs):
        return super().dispatch(*args, **kwargs)


class MicroscopeList(ListView):
    def get_queryset(self):
        # get all collections for current user and all other non-private collections
        qs = Microscope.objects.all()
        if self.request.user.is_authenticated:
            qs = qs | Microscope.objects.filter(owner=self.request.user)
        if 'owner' in self.kwargs:
            qs = qs.filter(owner__username=self.kwargs['owner'])
        return qs.order_by('-created')

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        if 'owner' in self.kwargs:
            if not get_user_model().objects.filter(username=self.kwargs['owner']).exists():
                raise Http404()
            context['owner'] = self.kwargs['owner']
        return context


class MicroscopeDeleteView(DeleteView):
    model = Microscope
    success_url = reverse_lazy('proteins:newmicroscope')

    def get_success_url(self):
        redirect_url = reverse_lazy('proteins:microscopes', kwargs={'owner': self.request.user})
        try:
            # check that this is an internal redirection
            resolve(redirect_url)
        except Exception:
            redirect_url = None
        return redirect_url or super().get_success_url()
