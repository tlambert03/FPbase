
import json
from django.views.generic import DetailView, CreateView
from django.http import HttpResponseRedirect
from django.core.mail import mail_admins
from ..models import Microscope, Camera, Light, Spectrum
from ..forms import MicroscopeForm
from .mixins import OwnableObject


class MicroscopeCreateView(OwnableObject, CreateView):
    ''' renders html for microscope creation page '''
    model = Microscope
    form_class = MicroscopeForm

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


class MicroscopeDetailView(DetailView):
    ''' renders html for microscope detail/spectrum page '''
    queryset = Microscope.objects.all().prefetch_related(
        'configs__filterplacement_set__filter',
        'configs__filters')

    def get_context_data(self, **kwargs):
        data = super().get_context_data(**kwargs)
        if not self.object.cameras.count():
            data['cameras'] = Camera.objects.all()
        if not self.object.lights.count():
            data['lights'] = Light.objects.all()
        data['probeslugs'] = Spectrum.objects.fluorlist()
        data['scopespectra'] = json.dumps(self.object.spectra_d3())
        # data['configs'] = json.dumps([oc.to_json() for oc in self.object.configs.all()])
        return data
