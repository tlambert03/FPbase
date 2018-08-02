import reversion
from django.views.generic import DetailView, CreateView, UpdateView
from django.views.generic.base import TemplateView
from django.views.decorators.cache import cache_page
from django.shortcuts import render, get_object_or_404
from django.http import HttpResponseRedirect, HttpResponse, JsonResponse
from django.forms.models import modelformset_factory
from django.contrib.auth.decorators import login_required
from django.contrib.admin.views.decorators import staff_member_required
from django.db import transaction
from django.db.models import Case, When
from django.db.models.functions import Substr
from django.core.mail import mail_managers
from ..models import Protein, State, Organism, BleachMeasurement
from ..forms import (ProteinForm, StateFormSet, StateTransitionFormSet,
                     BleachMeasurementForm)
from proteins.util.spectra import spectra2csv
from references.models import Reference  # breaks application modularity
from reversion.views import _RollBackRevisionView
from reversion.models import Version
import json


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

        similar = Protein.visible.filter(name__iexact='m' + self.object.name)
        similar = similar | Protein.visible.filter(name__iexact='monomeric' + self.object.name)
        similar = similar | Protein.visible.filter(name__iexact=self.object.name.lstrip('m'))
        similar = similar | Protein.visible.filter(name__iexact=self.object.name.lstrip('monomeric'))
        similar = similar | Protein.visible.filter(name__iexact=self.object.name.lower().lstrip('td'))
        similar = similar | Protein.visible.filter(name__iexact='td' + self.object.name)
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
                        sbj = 'Protein {} {}'.format(
                            self.object, self.get_form_type() + 'ed')
                        msg = "User: {}\nProtein: {}\n{}".format(
                            self.request.user.username,
                            self.object,
                            self.request.build_absolute_uri(
                                self.object.get_absolute_url()))
                        mail_managers(sbj, msg, fail_silently=True)
                    else:
                        self.object.status = 'approved'

                    self.object.save()
                    reversion.set_user(self.request.user)
                    reversion.set_comment('{} {} form'.format(
                        self.object, self.get_form_type()))
            else:
                context.update({
                    'states': states,
                })
                return self.render_to_response(context)

        return HttpResponseRedirect(self.get_success_url())


# https://medium.com/@adandan01/django-inline-formsets-example-mybook-420cc4b6225d
class ProteinCreateView(ProteinCreateUpdateMixin, CreateView):
    ''' renders html for protein submission page '''
    model = Protein
    form_class = ProteinForm

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
        form.instance.created_by = self.request.user
        return super().form_valid(form)


class ProteinUpdateView(ProteinCreateUpdateMixin, UpdateView):
    ''' renders html for protein submission page '''
    model = Protein
    form_class = ProteinForm

    def get_context_data(self, **kwargs):
        data = super().get_context_data(**kwargs)
        if self.request.POST:
            data['states'] = StateFormSet(
                self.request.POST, instance=self.object)
            data['states'].full_clean()  # why is this here?
        else:
            data['states'] = StateFormSet(instance=self.object)
            if self.object.primary_reference:
                data['form'].fields['reference_doi'].initial = self.object.primary_reference.doi
        return data

    def form_valid(self, form):
        form.instance.updated_by = self.request.user
        return super().form_valid(form)


@cache_page(60 * 10)
def protein_table(request):
    ''' renders html for protein table page  '''
    return render(
        request,
        'table.html',
        {
            "proteins": Protein.objects.all().prefetch_related(
                'states', 'states__bleach_measurements'),
            'request': request
        })


class ComparisonView(TemplateView):

    template_name = "compare.html"

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        if 'proteins' in self.kwargs:
            # the page was requested with slugs in the URL, maintain that order
            ids = kwargs.get('proteins', '').split(',')
            p = Case(*[When(slug=slug, then=pos) for pos, slug in enumerate(ids)])
            prots = Protein.objects.filter(slug__in=ids).prefetch_related('states__spectra').order_by(p)
        else:
            # try to use chronological order
            ids = self.request.session.get('comparison', [])
            prots = Protein.objects.filter(slug__in=ids).prefetch_related('states__spectra').order_by('primary_reference__year')
        if prots.exclude(seq__isnull=True).count() > 2:
            a, _ = prots.exclude(seq__isnull=True).to_tree('html')
            context['alignment'] = "\n".join(a.splitlines()[3:]).replace("FFEEE0", "FFFFFF")
        elif prots.count() == 2:
            if prots[0].seq and prots[1].seq:
                a = prots[0].seq.align_to(prots[1].seq)
                out = []
                for i, row in enumerate(str(a).splitlines()):
                    head = ''
                    if i % 4 == 0:
                        head = prots[0].name
                    elif i % 2 == 0:
                        head = prots[1].name
                    out.append("{:<16}{}".format(head, row))
                out.append('\n')
                context['alignment'] = "\n".join(out)
                context['mutations'] = a.mutations
        else:
            context['alignment'] = None
        context['proteins'] = prots
        refs = Reference.objects.filter(proteins__in=prots).distinct('id').order_by('id')
        context['references'] = sorted([r for r in refs], key=lambda x: x.year)
        spectra = []
        for prot in prots:
            for state in prot.states.all():
                spectra.extend(state.d3_dicts())
        context['spectra'] = json.dumps(spectra) if spectra else None
        return context


def protein_tree(request, organism):
    ''' renders html for protein table page  '''
    _, tree = Protein.objects.filter(parent_organism=organism).to_tree()
    return render(
        request,
        'tree.html',
        {
            "tree": tree.replace('\n', ''),
            'request': request
        })


def sequence_problems(request):
    ''' renders html for protein table page  '''
    mprobs = Protein.objects.annotate(sub=Substr('seq', 206, 1)).filter(sub='A', name__istartswith='m')
    mprobs = mprobs | Protein.objects.annotate(sub=Substr('seq', 207, 1)).filter(sub='A', name__istartswith='m')

    return render(
        request,
        'seq_problems.html',
        {
            "histags": Protein.objects.filter(seq__icontains='HHHHHH'),
            "noseqs": Protein.objects.filter(seq__isnull=True),
            "mprobs": mprobs,
            "nomet": Protein.objects.exclude(seq__isnull=True).exclude(seq__istartswith='M'),
            "noparent": Protein.objects.filter(parent_organism__isnull=True),
            'request': request
        })


def add_reference(request, slug=None):
    try:
        with reversion.create_revision():
            doi = request.POST.get('reference_doi')
            P = Protein.objects.get(slug=slug)
            ref, created = Reference.objects.get_or_create(doi=doi)
            P.references.add(ref)
            if not request.user.is_staff:
                P.status = 'pending'
                msg = "User: {}\nProtein: {}\nReference: {}, {}\n{}".format(
                    request.user.username, P, ref, ref.title,
                    request.build_absolute_uri(P.get_absolute_url()))
                mail_managers('Reference Added', msg, fail_silently=True)
            P.save()
            reversion.set_user(request.user)
            reversion.set_comment('Ref: {} added to {}'.format(ref, P))
        return JsonResponse({})
    except Exception:
        pass


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
                        mail_managers(
                            'Transition updated',
                            "User: {}\nProtein: {}\nForm: {}"
                            .format(request.user.username, obj, formset),
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
    BleachMeasurementFormSet = modelformset_factory(
        BleachMeasurement, BleachMeasurementForm, extra=1, can_delete=True)
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
                        mail_managers(
                            'BleachMeasurement Added',
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
    queryset = Organism.objects.all().prefetch_related('proteins__states')


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


