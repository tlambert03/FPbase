from django.views.generic import DetailView, ListView, CreateView, UpdateView
from django.shortcuts import render, redirect
from django.http import HttpResponseRedirect, HttpResponse, HttpResponseNotAllowed, HttpResponseBadRequest
from django.contrib import messages
from django.utils.text import slugify
from django.contrib.auth.decorators import login_required
from django.contrib.admin.views.decorators import staff_member_required
from django.core import serializers
from django.db import transaction
from django.http import JsonResponse
from django.contrib.postgres.search import TrigramSimilarity

from .models import Protein, State, ProteinCollection
from .forms import ProteinSubmitForm, ProteinUpdateForm, StateFormSet, StateTransitionFormSet

from references.models import Reference  # breaks application modularity # FIXME
import reversion
from reversion.views import _RollBackRevisionView
from reversion.models import Version


class CollectionList(ListView):
    def get_queryset(self):
        # get all collections for current user and all other non-private collections
        qs = ProteinCollection.objects.exclude(private=True)
        if self.request.user.is_authenticated:
            qs = qs | ProteinCollection.objects.filter(owner=self.request.user)
        if 'owner' in self.kwargs:
            qs = qs.filter(owner__username=self.kwargs['owner'])
        return qs.order_by('modified')

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        if 'owner' in self.kwargs:
            context['owner'] = self.kwargs['owner']
        return context


class CollectionDetail(DetailView):
    queryset = ProteinCollection.objects.all().prefetch_related('proteins', 'proteins__states', 'proteins__default_state')

    def get_context_data(self, **kwargs):
        # Call the base implementation first to get a context
        context = super().get_context_data(**kwargs)
        context['isowner'] = self.request.user == self.object.owner
        return context


class ProteinChartList(ListView):
    ''' renders html for interactive chart  '''
    template_name = 'ichart.html'
    model = Protein
    queryset = Protein.objects.filter(switch_type=Protein.BASIC).select_related('default_state')

    def get_context_data(self, **kwargs):
        # Call the base implementation first to get a context
        context = super().get_context_data(**kwargs)
        # Add in a QuerySet of all the books
        filtered_data = []
        filtered_data = serializers.serialize('json', self.get_queryset,
            fields=('name', 'default_state__em_max')
        )
        context['data'] = filtered_data
        return context


class ProteinDetailView(DetailView):
    ''' renders html for single protein page  '''
    queryset = Protein.objects.all().prefetch_related('states')

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
            try:
                version = Version.objects.get(id=kwargs['ver'])
                if int(version.object_id) == self.get_object().id:
                    return self.version_view(request, version, *args, **kwargs)
            except Version.DoesNotExist:
                pass
        return super().get(request, *args, **kwargs)

    def get_context_data(self, **kwargs):
        data = super().get_context_data(**kwargs)
        if not self.object.status == 'approved':
            data['last_approved'] = self.object.last_approved_version()
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
    form_class = ProteinSubmitForm
    # success_url --> used for redirect on success... by default shows new protein

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
    form_class = ProteinUpdateForm
    # success_url --> used for redirect on success... by default shows new protein

    def get_context_data(self, **kwargs):
        data = super().get_context_data(**kwargs)
        if self.request.POST:
            data['states'] = StateFormSet(self.request.POST, instance=self.object)
            data['states'].full_clean()
        else:
            data['states'] = StateFormSet(instance=self.object)
            data['states'].extra = 0
            if self.object.primary_reference:
                data['form'].fields['reference_doi'].initial = self.object.primary_reference.doi
        return data

    def form_valid(self, form):
        form.instance.updated_by = self.request.user  # login_required in url.py
        return super().form_valid(form)


class TransitionUpdateView(UpdateView):
    model = Protein

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        if self.request.POST:
            context['transition_form'] = StateTransitionFormSet(self.request.POST,  instance=self.object)
        else:
            context['transition_form'] = StateTransitionFormSet()
        return context

    def form_valid(self, form):
        context = self.get_context_data()
        transition_form = context['transition_form']
        if transition_form.is_valid():
            self.object = form.save()
            transition_form.instance = self.object
            transition_form.save()
            return HttpResponseRedirect(self.object.get_absolute_url())
        else:
            return self.render_to_response(self.get_context_data(form=form))


def protein_table(request):
    ''' renders html for protein table page  '''
    #  return render(request, 'table.html', {"states": State.objects.notdark().select_related('protein').prefetch_related('bleach_measurements')})
    return render(request, 'table.html', {"proteins": Protein.objects.all().prefetch_related('states', 'states__bleach_measurements')})


def protein_search(request):
    ''' renders html for protein search page  '''
    from .filters import ProteinFilter

    if request.GET:
        f = ProteinFilter(request.GET, queryset=Protein.objects.select_related('default_state').prefetch_related('states').order_by('default_state__em_max'))

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
            P.save()
            reversion.set_user(request.user)
            reversion.set_comment('Ref: {} added to {}'.format(ref, P))
        return JsonResponse({})
    except Exception:
        pass


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
    name = request.POST.get('name', None)
    try:
        prot = Protein.objects.get(slug=slugify(name.replace(' ', '')))
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
def collection_remove(request):

    if not request.is_ajax():
        return HttpResponseNotAllowed([])

    try:
        protein = int(request.POST["target_protein"])
        collection = int(request.POST["target_collection"])
    except (KeyError, ValueError):
        return HttpResponseBadRequest()

    col = ProteinCollection.objects.get(id=collection)
    if not col.owner == request.user:
        return HttpResponseNotAllowed([])

    col.proteins.remove(protein)

    response = {
        'status': 'deleted',
    }

    return JsonResponse(response)

