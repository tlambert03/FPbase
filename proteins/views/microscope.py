import json
from django.views.generic import DetailView, CreateView, DeleteView, ListView, UpdateView
from django.http import HttpResponseRedirect, Http404
from django.core.mail import mail_admins
from django.urls import reverse_lazy
from django.contrib.auth import get_user_model
from django.utils.decorators import method_decorator
from django.views.decorators.clickjacking import xframe_options_exempt

from ..models import Microscope, Camera, Light, Spectrum, OpticalConfig
from ..forms import MicroscopeForm, OpticalConfigFormSet
from .mixins import OwnableObject


class MicroscopeCreateUpdateMixin:

    def get_form_type(self):
        return self.request.resolver_match.url_name

    def form_valid(self, form):
        self.attach_owner(form)
        if not self.request.user.is_staff:
            mail_admins('Microscope Created',
                        "User: {}\nCollection: {}\n{}".format(
                            self.request.user.username,
                            self.object,
                            self.request.build_absolute_uri(
                                self.object.get_absolute_url())),
                        fail_silently=True)
        return HttpResponseRedirect(self.get_success_url())


class MicroscopeCreateView(MicroscopeCreateUpdateMixin, OwnableObject, CreateView):
    ''' renders html for microscope creation page '''
    model = Microscope
    form_class = MicroscopeForm

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


class MicroscopeUpdateView(MicroscopeCreateUpdateMixin, OwnableObject, UpdateView):
    model = Microscope
    form_class = MicroscopeForm

    def form_valid(self, form):
        form.instance.updated_by = self.request.user
        return super().form_valid(form)

    def get_context_data(self, **kwargs):
        data = super().get_context_data(**kwargs)
        if self.request.POST:
            data['optical_configs'] = OpticalConfigFormSet(self.request.POST,
                                                           instance=self.object)
        else:
            inst = self.object
            qs = inst.optical_configs.prefetch_related(
                'ex_filters', 'em_filters', 'bs_filters', 'filters')
            formdata = {
                'form-TOTAL_FORMS': qs.count() + 1,
                'form-INITIAL_FORMS': 0,
            }
            for i, oc in enumerate(qs):
                formdata['form-{}-ex_filters'.format(i)] = [f.id for f in oc.ex_filters.all()]
                formdata['form-{}-em_filters'.format(i)] = [f.id for f in oc.em_filters.all()]
                formdata['form-{}-bs_filters'.format(i)] = [f.id for f in oc.bs_filters.all()]
                formdata['form-{}-name'.format(i)] = getattr(oc, 'name', None)
                formdata['form-{}-laser'.format(i)] = getattr(oc, 'laser', None)
                formdata['form-{}-light'.format(i)] = getattr(oc.light, 'id', None)
                formdata['form-{}-camera'.format(i)] = getattr(oc.camera, 'id', None)
                formdata['form-{}-invert_bs'.format(i)] = getattr(oc, 'inverted_bs', False)
            form = OpticalConfigFormSet(formdata, instance=inst, prefix='form',
                                        queryset=OpticalConfig.objects.all())
            data['optical_configs'] = form
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
        data['probeslugs'] = Spectrum.objects.fluorlist()
        data['scopespectra'] = json.dumps(self.object.spectra_d3())
        # data['optical_configs'] = json.dumps([oc.to_json() for oc in self.object.optical_configs.all()])
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

    # def get_success_url(self):
    #     redirect_url = reverse_lazy('proteins:collections', kwargs={'owner': self.request.user})
    #     try:
    #         # check that this is an internal redirection
    #         resolve(redirect_url)
    #     except Exception:
    #         redirect_url = None
    #     return redirect_url or super().get_success_url()

