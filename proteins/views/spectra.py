import json
from django.conf import settings
from django.core.mail import EmailMessage
from django.core.cache import cache
from django.http import Http404, HttpResponse, JsonResponse
from django.urls import reverse
from django.views.generic import CreateView
from django.views.decorators.cache import cache_page
from django.views.decorators.vary import vary_on_cookie
from django.shortcuts import render
from django.template.defaultfilters import slugify
from ..models import Spectrum, Filter, State, Protein
from ..forms import SpectrumForm
from ..util.spectra import spectra2csv
from ..util.importers import add_filter_to_database
from ..util.helpers import forster_list


@cache_page(60 * 10)
@vary_on_cookie
def protein_spectra(request, slug=None):
    ''' renders html for protein spectra page  '''
    template = 'spectra.html'

    if request.is_ajax() and slug is not None:
        ''' slug represents the slug of the OWNER of the spectra...
        for instance, a protein state.
        function returns an array containing ALL of the spectra
        belonging to that owner
        '''
        owner = Spectrum.objects.filter_owner(slug)
        if owner.count():
            return JsonResponse({'spectra': [s.d3dict() for s in owner]})
        else:
            raise Http404

    if request.method == 'GET':
        # spectra = Spectrum.objects.owner_slugs()
        spectra = Spectrum.objects.sluglist()

        return render(request, template, {
            'spectra_options': json.dumps(spectra)}
        )


def fret_chart(request):
    ''' renders html for protein spectra page  '''
    template = 'fret.html'

    if request.method == 'GET':

        if request.is_ajax():
            L = cache.get('forster_list')
            if not L:
                L = forster_list()
                cache.set('forster_list', L, 60 * 60 * 24 * 7)
            return JsonResponse({'data': L})

        else:
            slugs = State.objects.filter(spectra__subtype='ex').values('slug', 'protein__name', 'spectra__category', 'spectra__subtype')
            slugs = [{'slug': x['slug'], 'category': x['spectra__category'], 'subtype': x['spectra__subtype'], 'name': x['protein__name']} for x in slugs]
            return render(request, template, {
                'probeslugs': slugs}
            )


class SpectrumCreateView(CreateView):
    model = Spectrum
    form_class = SpectrumForm

    def get_success_url(self, **kwargs):
        if self.object.category == Spectrum.PROTEIN:
            return self.object.owner.get_absolute_url()
        else:
            return "{}?s={}".format(reverse('proteins:spectra'), self.object.owner.slug)

    def get_initial(self):
        init = super().get_initial()
        if self.kwargs.get('slug', False):
            try:
                p = Protein.objects.get(slug=self.kwargs.get('slug'))
                init['owner_state'] = p.default_state
                init['category'] = Spectrum.PROTEIN
            except Exception:
                pass
        return init

    def get_form_kwargs(self):
        kwargs = super().get_form_kwargs()
        kwargs['user'] = self.request.user
        return kwargs

    def get_form(self, form_class=None):
        form = super().get_form()

        if self.kwargs.get('slug', False):
            try:
                from django import forms
                p = Protein.objects.get(slug=self.kwargs.get('slug'))
                form.fields['owner_state'] = forms.ModelChoiceField(
                    required=True, label='Protein', empty_label=None,
                    queryset=State.objects.filter(protein=p).select_related('protein'))
            except Exception:
                pass
        return form

    def form_valid(self, form):
        # This method is called when valid form data has been POSTed.
        # It should return an HttpResponse.
        i = super().form_valid(form)
        if not self.request.user.is_staff:
            EmailMessage(
                '[FPbase] Spectrum submitted: %s' % form.cleaned_data['owner'],
                self.request.build_absolute_uri(form.instance.get_admin_url()),
                to=[a[1] for a in settings.ADMINS],
                headers={'X-Mailgun-Track': 'no'}
            ).send()
        return i

    def get_context_data(self, **kwargs):
        # Call the base implementation first to get a context
        context = super().get_context_data(**kwargs)
        return context


def spectra_csv(request):
    try:
        idlist = [int(x) for x in request.GET.get('q', '').split(',') if x]
        spectralist = Spectrum.objects.filter(id__in=idlist)
        if spectralist:
            return spectra2csv(spectralist)
    except Exception:
        return HttpResponse('malformed spectra csv request')
    else:
        return HttpResponse('malformed spectra csv request')


def filter_import(request, brand):
    part = request.POST['part']
    newObjects = []
    errors = []
    response = {'status': 0}

    try:
        Filter.objects.get(slug=slugify(brand + ' ' + part))
        response['message'] = '%s is already in the database' % part
        return JsonResponse(response)
    except Filter.DoesNotExist:
        pass

    try:
        newObjects, errors = add_filter_to_database(brand, part, request.user)
    except Exception as e:
        response['message'] = str(e)

    if newObjects:
        spectrum = newObjects[0]
        response = {
            'status': 1,
            'objects': spectrum.name,
            'spectra_options': json.dumps({
                'category': spectrum.category,
                'subtype': spectrum.subtype,
                'slug': spectrum.owner.slug,
                'name': spectrum.owner.name,
            })
        }
    elif errors:
        try:
            if errors[0][1].as_data()['owner'][0].code == 'owner_exists':
                response['message'] = '%s already appears to be imported' % part
        except Exception:
            pass
    return JsonResponse(response)
