from django.views.generic import DetailView, ListView, CreateView, UpdateView, DeleteView
from django.shortcuts import render, redirect, get_object_or_404
from django.http import HttpResponseRedirect, HttpResponse, HttpResponseNotAllowed, HttpResponseBadRequest
from django.contrib import messages
from django.forms.models import modelformset_factory
from django.utils.text import slugify
from django.contrib.auth.decorators import login_required
from django.contrib.admin.views.decorators import staff_member_required
from django.db import transaction
from django.http import JsonResponse
from django.contrib.postgres.search import TrigramSimilarity
from django import forms
from django.urls import resolve, reverse_lazy
from django.core.mail import mail_managers, mail_admins

import json
from .models import Protein, State, ProteinCollection, Organism, BleachMeasurement
from .forms import (ProteinForm, StateFormSet, StateTransitionFormSet,
                    CollectionForm, BleachMeasurementForm)
from .filters import ProteinFilter

from references.models import Reference  # breaks application modularity # FIXME
import reversion
from reversion.views import _RollBackRevisionView
from reversion.models import Version


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
        data['similar'] = similar.exclude(id=self.object.id)
        return data


class ProteinCreateUpdateMixin:

    def get_form_type(self):
        return self.request.resolver_match.url_name

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


# class TransitionUpdateView(UpdateView):
#     model = Protein

#     def get_context_data(self, **kwargs):
#         context = super().get_context_data(**kwargs)
#         if self.request.POST:
#             context['transition_form'] = StateTransitionFormSet(self.request.POST,  instance=self.object)
#         else:
#             context['transition_form'] = StateTransitionFormSet()
#         return context

#     def form_valid(self, form):
#         context = self.get_context_data()
#         transition_form = context['transition_form']
#         if transition_form.is_valid():
#             self.object = form.save()
#             transition_form.instance = self.object
#             transition_form.save()
#             return HttpResponseRedirect(self.object.get_absolute_url())
#         else:
#             return self.render_to_response(self.get_context_data(form=form))


def protein_table(request):
    ''' renders html for protein table page  '''
    return render(request, 'table.html', {"proteins": Protein.objects.all().prefetch_related('states', 'states__bleach_measurements')})


def protein_spectra(request, slug=None):
    ''' renders html for protein spectra page  '''
    if request.is_ajax() and slug is not None:
        protein = Protein.objects.get(slug=slug)
        return JsonResponse({'spectra': protein.spectra_json()})
    if request.method == 'GET':
        qs = Protein.objects.with_spectra().values_list('slug', 'name')
        widget = forms.Select(attrs={'class': 'form-control custom-select', 'id': 'ProteinSelect'})
        choicefield = forms.ChoiceField(choices=qs, widget=widget)

    return render(request, 'spectra.html', {'widget': choicefield.widget.render('ProteinSelect', '')})


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


@login_required
def add_organism(request):
    if not request.is_ajax():
        return HttpResponseNotAllowed([])
    try:
        tax_id = request.POST.get('taxonomy_id', None)
        if not tax_id:
            raise Exception()
        org, created = Organism.objects.get_or_create(id=tax_id)
        if created:
            if request.user.is_authenticated:
                org.created_by = request.user
                org.save()
            if not request.user.is_staff:
                mail_managers('Organism Added',
                    "User: {}\nOrganism: {}\n{}".format(
                        request.user.username, org,
                        request.build_absolute_uri(org.get_absolute_url())),
                    fail_silently=True)
        return JsonResponse({'status': 'success', 'created': created,
                             'org_id': org.id, 'org_name': org.scientific_name})
    except Exception:
        return JsonResponse({'status': 'failed'})


@staff_member_required
def approve_protein(request, slug=None):
    try:
        P = Protein.objects.get(slug=slug)
        if not P.status == 'pending':
            return JsonResponse({})

        # get rid of previous unapproved version
        try:
            if P.versions.first().field_dict['status'] == 'pending':
                P.versions.first().delete()
        except Exception:
            pass

        with reversion.create_revision():
            P.status = 'approved'
            P.save()

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


def validate_proteinname(request):
    if not request.is_ajax():
        return HttpResponseNotAllowed([])

    name = request.POST.get('name', None)
    slug = request.POST.get('slug', None)
    try:
        prot = Protein.objects.get(slug=slugify(name.replace(' ', '').replace('monomeric', 'm')))
        if slug and prot.slug == slug:
            data = {
                'is_taken': False,
            }
        else:
            data = {
                'is_taken': True,
                'id': prot.id,
                'url': prot.get_absolute_url(),
                'name': prot.name,
            }
    except Protein.DoesNotExist:
        data = {
            'is_taken': False,
        }
    return JsonResponse(data)


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


class CollectionList(ListView):
    def get_queryset(self):
        # get all collections for current user and all other non-private collections
        qs = ProteinCollection.objects.exclude(private=True)
        if self.request.user.is_authenticated:
            qs = qs | ProteinCollection.objects.filter(owner=self.request.user)
        if 'owner' in self.kwargs:
            qs = qs.filter(owner__username=self.kwargs['owner'])
        return qs.order_by('-created')

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        if 'owner' in self.kwargs:
            context['owner'] = self.kwargs['owner']
        return context


def serialized_proteins_response(queryset, format='json', filename='FPbase_proteins'):
    from proteins.api.serializers import ProteinSerializer as PS
    PS.Meta.on_demand_fields = ()
    serializer = PS(queryset, many=True)
    if format == 'json':
        from rest_framework.renderers import JSONRenderer as rend
        response = JsonResponse(serializer.data, safe=False)
    elif format == 'csv':
        from rest_framework_csv.renderers import CSVStreamingRenderer as rend
        from django.http import StreamingHttpResponse
        response = StreamingHttpResponse(rend().render(serializer.data), content_type='text/csv')
        response['Content-Disposition'] = 'attachment; filename="{}.csv"'.format(filename)
    return response


class CollectionDetail(DetailView):
    queryset = ProteinCollection.objects.all().prefetch_related('proteins', 'proteins__states', 'proteins__default_state')

    def get(self, request, *args, **kwargs):
        format = request.GET.get('format', '').lower()
        if format in ('json', 'csv'):
            col = self.get_object()
            return serialized_proteins_response(col.proteins.all(), format,
                    filename=slugify(col.name))
        return super().get(request, *args, **kwargs)

    def get_context_data(self, **kwargs):
        # Call the base implementation first to get a context
        context = super().get_context_data(**kwargs)
        context['isowner'] = self.request.user == self.object.owner
        return context

    def render_to_response(self, *args, **kwargs):
        if not self.request.user.is_superuser:
            if self.object.private and (self.object.owner != self.request.user):
                return render(self.request, 'proteins/private_collection.html', {'foo': 'bar'})
        return super().render_to_response(*args, **kwargs)


@login_required
def collection_remove(request):
    if not request.is_ajax():
        return HttpResponseNotAllowed([])
    try:
        protein = int(request.POST["target_protein"])
        collection = int(request.POST["target_collection"])
    except (KeyError, ValueError):
        return HttpResponseBadRequest()

    col = get_object_or_404(ProteinCollection, id=collection)

    if not col.owner == request.user:
        return HttpResponseNotAllowed([])
    col.proteins.remove(protein)
    response = {
        'status': 'deleted',
    }
    return JsonResponse(response)


@login_required
def add_to_collection(request):
    if not request.is_ajax():
        return HttpResponseNotAllowed([])

    if request.method == 'GET':
        qs = ProteinCollection.objects.filter(owner=request.user)
        widget = forms.Select(attrs={'class': 'form-control custom-select', 'id': 'collectionSelect'})
        choicefield = forms.ChoiceField(choices=qs.values_list('id', 'name'), widget=widget)

        members = []
        if request.GET.get('id'):
            try:
                qs = qs.filter(proteins=int(request.GET.get('id')))
                members = [(item.name, item.get_absolute_url()) for item in qs]
            except Exception as e:
                print(e)
                pass

        response = {
            'widget': choicefield.widget.render('collectionChoice', ''),
            'members': json.dumps(members),
        }
        return JsonResponse(response)

    elif request.method == 'POST':
        try:
            collection = ProteinCollection.objects.get(id=request.POST.get('collectionChoice'))
            collection.proteins.add(int(request.POST.get('protein')))
            status = 'success'
        except Exception as e:
            status = 'error'
            print(e)
        return JsonResponse({'status': status})

    return HttpResponseNotAllowed([])


class CollectionCreateView(CreateView):
    model = ProteinCollection
    form_class = CollectionForm

    def get_form_kwargs(self):
        kwargs = super().get_form_kwargs()
        kwargs['user'] = self.request.user
        if len(self.request.POST.getlist('protein')):
            kwargs['proteins'] = self.request.POST.getlist('protein')
        elif self.request.POST.get('dupcollection', False):
            id = self.request.POST.get('dupcollection')
            kwargs['proteins'] = [p.id for p in ProteinCollection.objects.get(id=id).proteins.all()]
        return kwargs

    def form_valid(self, form):
        self.object = form.save(commit=False)
        self.object.owner = self.request.user
        self.object.save()
        if getattr(form, 'proteins', None):
            self.object.proteins.add(*form.proteins)
        if not self.request.user.is_staff:
            mail_admins('Collection Created',
                "User: {}\nCollection: {}\n{}".format(
                    self.request.user.username,
                    self.object,
                    self.request.build_absolute_uri(self.object.get_absolute_url())),
                fail_silently=True)
        return HttpResponseRedirect(self.get_success_url())

    def get_success_url(self):
        redirect_url = self.request.POST.get('next') or self.request.GET.get('next', None)
        try:
            # check that this is an internal redirection
            resolve(redirect_url)
        except Exception:
            redirect_url = None
        return redirect_url or super().get_success_url()

    def get_context_data(self, **kwargs):
        data = super().get_context_data(**kwargs)
        if self.request.POST.get('colname', False):
            data['colname'] = self.request.POST.get('colname')
        return data


class CollectionUpdateView(UpdateView):
    model = ProteinCollection
    form_class = CollectionForm

    def get_form_kwargs(self):
        kwargs = super().get_form_kwargs()
        kwargs['user'] = self.request.user
        return kwargs


class CollectionDeleteView(DeleteView):
    model = ProteinCollection
    success_url = reverse_lazy('proteins:collections')

    def get_success_url(self):
        redirect_url = reverse_lazy('proteins:collections', kwargs={'owner': self.request.user})
        try:
            # check that this is an internal redirection
            resolve(redirect_url)
        except Exception:
            redirect_url = None
        return redirect_url or super().get_success_url()
