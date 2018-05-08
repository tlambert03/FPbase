from django.views.generic import DetailView, CreateView, UpdateView
from django.shortcuts import render, redirect, get_object_or_404
from django.http import HttpResponseRedirect, HttpResponse, Http404
from django.forms.models import modelformset_factory
from django.contrib.auth.decorators import login_required
from django.contrib.admin.views.decorators import staff_member_required
from django.db import transaction
from django.http import JsonResponse
from django.contrib.postgres.search import TrigramSimilarity
from django.core.mail import mail_managers
from django.core.exceptions import ValidationError
import json

from ..models import Protein, State, Organism, BleachMeasurement, Spectrum
from ..forms import (ProteinForm, StateFormSet, StateTransitionFormSet,
                     BleachMeasurementForm, SpectrumForm)
from ..filters import ProteinFilter

from references.models import Reference  # breaks application modularity # FIXME
import reversion
from reversion.views import _RollBackRevisionView
from reversion.models import Version
from ..util.importers import import_chroma_spectra, import_semrock_spectra
from django.core.validators import URLValidator


class ProteinDetailView(DetailView):
    ''' renders html for single protein page  '''
    queryset = Protein.visible.all().prefetch_related('states')

    def version_view(self, request, version, *args, **kwargs):
        try:
            with transaction.atomic(using=version.db):
                # Revert the revision.
                version.revision.revert(delete=True)
                # Run the normal changeform view.
                self.object = self.get_object()
                context = self.get_context_data(object=self.object)
                context['version'] = version
                response = self.render_to_response(context)
                response.render()  # eager rendering of response is necessary before db rollback
                raise _RollBackRevisionView(response)
        except _RollBackRevisionView as ex:
            return ex.response

    def get(self, request, *args, **kwargs):
        if 'rev' in kwargs:
            try:
                rev = int(kwargs['rev'])  # has to be int or indexing will fail
            except Exception:
                rev = 0
            if rev > 0:
                versions = Version.objects.get_for_object(self.get_object())
                version = versions[min(versions.count() - 1, rev)]
                return self.version_view(request, version, *args, **kwargs)
        elif 'ver' in kwargs:
            version = get_object_or_404(Version, id=kwargs['ver'])
            if int(version.object_id) == self.get_object().id:
                return self.version_view(request, version, *args, **kwargs)
            # TODO:  ELSE WHAT??
        return super().get(request, *args, **kwargs)

    def get_context_data(self, **kwargs):
        data = super().get_context_data(**kwargs)
        if not self.object.status == 'approved':
            data['last_approved'] = self.object.last_approved_version()

        similar = Protein.visible.filter(name__iexact='m'+self.object.name)
        similar = similar | Protein.visible.filter(name__iexact='monomeric'+self.object.name)
        similar = similar | Protein.visible.filter(name__iexact=self.object.name.lstrip('m'))
        similar = similar | Protein.visible.filter(name__iexact=self.object.name.lstrip('monomeric'))
        similar = similar | Protein.visible.filter(name__iexact=self.object.name.lower().lstrip('td'))
        similar = similar | Protein.visible.filter(name__iexact='td'+self.object.name)
        data['similar'] = similar.exclude(id=self.object.id)
        return data


class ProteinCreateUpdateMixin:

    def get_form_type(self):
        return self.request.resolver_match.url_name

    # from django.contrib import messages
    # TODO: pull out the good message parts... remove the moderation parts
    # def moderate(self, obj):
    #     from moderation.helpers import automoderate
    #     status = automoderate(obj, self.request.user)
    #     if isinstance(obj, Protein):
    #         if status == 2:  # PENDING
    #             if self.get_form_type() == 'update':
    #                 messages.add_message(self.request, messages.INFO,
    #                     'Your update to {} has been submitted and will appear after moderation.'.format(obj))
    #             else:
    #                 messages.add_message(self.request, messages.INFO,
    #                     'Thank you for submitting {}.  It will appear shortly, after moderation.'.format(obj))
    #         elif status == 1:  # APPROVED
    #             if self.get_form_type() == 'update':
    #                 messages.add_message(self.request, messages.SUCCESS,
    #                     'Your update to {} has been approved.'.format(obj))
    #             else:
    #                 messages.add_message(self.request, messages.SUCCESS,
    #                     'Thank you for submitting {}!'.format(obj))
    #         else:  # REJECTED
    #             messages.add_message(self.request, messages.ERROR,
    #                 '{} rejected.  Please contact us with questions.'.format(self.get_form_type().title()))

    def form_valid(self, form):
        # This method is called when valid form data has been POSTed.
        context = self.get_context_data()
        states = context['states']

        with transaction.atomic():
            # only save the form if all the states are also valid
            if states.is_valid():
                with reversion.create_revision():
                    self.object = form.save()
                    doi = form.cleaned_data.get('reference_doi')
                    if doi:
                        ref, created = Reference.objects.get_or_create(doi=doi)
                        self.object.primary_reference = ref
                    else:
                        self.object.primary_reference = None

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

                    if not self.request.user.is_staff:
                        self.object.status = 'pending'
                        mail_managers('Protein {} {}'.format(self.object, self.get_form_type() + 'ed'),
                            "User: {}\nProtein: {}\n{}".format(
                                self.request.user.username,
                                self.object,
                                self.request.build_absolute_uri(self.object.get_absolute_url())),
                            fail_silently=True)
                    else:
                        self.object.status = 'approved'

                    self.object.save()
                    reversion.set_user(self.request.user)
                    reversion.set_comment('{} {} form'.format(self.object, self.get_form_type()))
            else:
                context.update({
                    'states': states,
                })
                return self.render_to_response(context)

        return HttpResponseRedirect(self.get_success_url())


# help from https://medium.com/@adandan01/django-inline-formsets-example-mybook-420cc4b6225d
class ProteinCreateView(ProteinCreateUpdateMixin, CreateView):
    ''' renders html for protein submission page '''
    model = Protein
    form_class = ProteinForm
    # success_url --> used for redirect on success... by default shows new protein

    def get_form(self, *args, **kwargs):
        form = super().get_form(*args, **kwargs)
        # make reference doi required when creating a protein
        form.fields['reference_doi'].required = True
        return form

    def get_context_data(self, **kwargs):
        data = super().get_context_data(**kwargs)
        if self.request.POST:
            data['states'] = StateFormSet(self.request.POST)
        else:
            data['states'] = StateFormSet()
        return data

    def form_valid(self, form):
        form.instance.created_by = self.request.user  # login_required in url.py
        return super().form_valid(form)


class ProteinUpdateView(ProteinCreateUpdateMixin, UpdateView):
    ''' renders html for protein submission page '''
    model = Protein
    form_class = ProteinForm
    # success_url --> used for redirect on success... by default shows new protein

    def get_context_data(self, **kwargs):
        data = super().get_context_data(**kwargs)
        if self.request.POST:
            data['states'] = StateFormSet(self.request.POST, instance=self.object)
            data['states'].full_clean()  # why is this here?
        else:
            data['states'] = StateFormSet(instance=self.object)
            if self.object.primary_reference:
                data['form'].fields['reference_doi'].initial = self.object.primary_reference.doi
        return data

    def form_valid(self, form):
        form.instance.updated_by = self.request.user  # login_required in url.py
        return super().form_valid(form)


def protein_table(request):
    ''' renders html for protein table page  '''
    return render(request, 'table.html', {"proteins": Protein.objects.all().prefetch_related('states', 'states__bleach_measurements')})


def protein_spectra(request, slug=None):
    ''' renders html for protein spectra page  '''
    template = 'spectra2.html'

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


class SpectrumCreateView(CreateView):
    model = Spectrum
    form_class = SpectrumForm


def protein_search(request):
    ''' renders html for protein search page  '''

    if request.GET:
        f = ProteinFilter(request.GET, queryset=Protein.objects.select_related('default_state')
                        .prefetch_related('states').order_by('default_state__em_max'))

        # if no hits, but name was provided... try trigram search
        if len(f.qs) == 0:
            name = None
            if 'name__icontains' in f.form.data:
                name = f.form.data['name__icontains']
            elif 'name__iexact' in f.form.data:
                name = f.form.data['name__iexact']
            if name:
                f.recs = Protein.objects.annotate(similarity=TrigramSimilarity('name', name)).filter(similarity__gt=0.2).order_by('-similarity')
        if len(f.qs) == 1:
            return redirect(f.qs.first())
    else:
        f = ProteinFilter(request.GET, queryset=Protein.objects.none())
    return render(request, 'proteins/protein_filter.html', {'filter': f})


def add_reference(request, slug=None):
    try:
        with reversion.create_revision():
            doi = request.POST.get('reference_doi')
            P = Protein.objects.get(slug=slug)
            ref, created = Reference.objects.get_or_create(doi=doi)
            P.references.add(ref)
            if not request.user.is_staff:
                P.status = 'pending'
                mail_managers('Reference Added',
                    "User: {}\nProtein: {}\nReference: {}, {}\n{}".format(
                        request.user.username, P, ref, ref.title,
                        request.build_absolute_uri(P.get_absolute_url())),
                    fail_silently=True)
            P.save()
            reversion.set_user(request.user)
            reversion.set_comment('Ref: {} added to {}'.format(ref, P))
        return JsonResponse({})
    except Exception:
        pass


def filter_import(request, brand):
    if brand.lower() == 'chroma':
        importer = import_chroma_spectra
    elif brand.lower() == 'semrock':
        importer = import_semrock_spectra
    else:
        raise ValueError('unknown brand')

    part = request.POST['part']
    newObjects = []
    errors = []
    response = {'status': 0}
    try:
        urlv = URLValidator()
        urlv(part)
        try:
            newObjects, errors = importer(url=part)
        except Exception as e:
            response['message'] = str(e)
    except ValidationError:
        try:
            newObjects, errors = importer(part=part)
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


@staff_member_required
def revert_version(request, ver=None):
    try:
        version = Version.objects.get(id=ver)
        version.revision.revert(delete=True)
        return JsonResponse({})
    except Exception:
        pass


@login_required
def update_transitions(request, slug=None):
    template_name = 'proteins/forms/_transition_form.html'
    obj = Protein.objects.get(slug=slug)
    if request.method == 'POST':
        formset = StateTransitionFormSet(request.POST, instance=obj)
        if formset.is_valid():
            with transaction.atomic():
                with reversion.create_revision():
                    formset.save()
                    if request.user.is_staff:
                        obj.status = 'approved'
                    else:
                        obj.status = 'pending'
                        mail_managers('Transition updated',
                            "User: {}\nProtein: {}\nForm: {}".format(request.user.username, obj, formset),
                            fail_silently=True)
                    obj.save()
                    reversion.set_user(request.user)
                    reversion.set_comment('Transitions edited on {}'.format(obj))
            return HttpResponse(status=200)
        else:
            response = render(request, template_name, {'transition_form': formset}, status=422)
            return response
    else:
        formset = StateTransitionFormSet(instance=obj)
        formset.form.base_fields['from_state'].queryset = State.objects.filter(protein=obj)
        formset.form.base_fields['to_state'].queryset = State.objects.filter(protein=obj)
        return render(request, template_name, {'transition_form': formset})


@login_required
def protein_bleach_formsets(request, slug):
    template_name = 'proteins/protein_bleach_form.html'
    BleachMeasurementFormSet = modelformset_factory(BleachMeasurement,
                        BleachMeasurementForm, extra=1, can_delete=True)
    protein = get_object_or_404(Protein, slug=slug)
    qs = BleachMeasurement.objects.filter(state__protein=protein)
    if request.method == 'POST':
        formset = BleachMeasurementFormSet(request.POST, queryset=qs)
        formset.form.base_fields['state'].queryset = State.objects.filter(protein__slug=slug)
        if formset.is_valid():
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

                    if not request.user.is_staff:
                        protein.status = 'pending'
                        mail_managers('BleachMeasurement Added',
                            "User: {}\nProtein: {}\n{}".format(
                                request.user.username, protein,
                                request.build_absolute_uri(protein.get_absolute_url())),
                            fail_silently=True)
                    else:
                        protein.status = 'approved'

                    protein.save()
                    reversion.set_user(request.user)
                    reversion.set_comment('Updated bleach measurement on {}'.format(protein))
            return HttpResponseRedirect(protein.get_absolute_url())
        else:
            return render(request, template_name, {'formset': formset, 'protein': protein})
    else:
        formset = BleachMeasurementFormSet(queryset=qs)
        formset.form.base_fields['state'].queryset = State.objects.filter(protein__slug=slug)
    return render(request, template_name, {'formset': formset, 'protein': protein})


class OrganismDetailView(DetailView):
    ''' renders html for single reference page  '''
    queryset = Organism.objects.all().prefetch_related('proteins')
