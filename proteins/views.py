from django.views.generic import DetailView, ListView, CreateView, UpdateView
from django.shortcuts import render, redirect
from django.http import HttpResponseRedirect
from django.contrib import messages
from django.core import serializers
from django.db import transaction
from django.http import JsonResponse
from django.contrib.postgres.search import TrigramSimilarity

from .models import Protein, State
from .forms import ProteinSubmitForm, ProteinUpdateForm, StateFormSet, StateUpdateFormSet

from references.models import Reference  # breaks application modularity # FIXME

from moderation import constants

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

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        try:
            if self.moderated_object.status == constants.MODERATION_STATUS_PENDING:
                context['pending_protein'] = self.object.moderated_object.changed_object
        except Exception:
            pass
        try:
            context['pending_states'] = [s.moderated_object.changed_object for s in State.unmoderated_objects.filter(protein=self.object).all()]

        except Exception:
            pass
        return context


class ProteinCreateUpdateMixin:

    def get_form_type(self):
        return self.request.resolver_match.url_name

    def moderate(self, obj):
        from moderation.helpers import automoderate
        status = automoderate(obj, self.request.user)
        if isinstance(obj, Protein):
            if status == 2:  # PENDING
                if self.get_form_type() == 'update':
                    messages.add_message(self.request, messages.INFO,
                        'Your update to {} has been submitted and will appear after moderation.'.format(obj))
                else:
                    messages.add_message(self.request, messages.INFO,
                        'Thank you for submitting {}.  It will appear shortly, after moderation.'.format(obj))
            elif status == 1:  # APPROVED
                if self.get_form_type() == 'update':
                    messages.add_message(self.request, messages.SUCCESS,
                        'Your update to {} has been approved.'.format(obj))
                else:
                    messages.add_message(self.request, messages.SUCCESS,
                        'Thank you for submitting {}!'.format(obj))
            else:  # REJECTED
                messages.add_message(self.request, messages.ERROR,
                    '{} rejected.  Please contact us with questions.'.format(self.get_form_type().title()))

    def form_valid(self, form):
        # This method is called when valid form data has been POSTed.
        context = self.get_context_data()
        states = context['states']

        with transaction.atomic():
            # only save the form if all the states are also valid
            if states.is_valid():
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
                    self.moderate(s)
                for s in states.deleted_objects:
                    s.delete()
                    # FIXME: state deletions currently unmoderated?
                    # self.moderate(s)

                self.object.save()
                self.moderate(self.object)
            else:
                context.update({
                    'states': states
                })
                return self.render_to_response(context)

        # if this is an update form, just pass on the sucess_url
        if self.get_form_type() == 'update':
            return HttpResponseRedirect(self.get_success_url())
        # otherwise we need to know whether there is a moderated protein to show or now
        else:
            if self.object in Protein.objects.all():
                return HttpResponseRedirect(self.get_success_url())
            # if not: just go back to the submit page
            else:
                return HttpResponseRedirect(self.request.path_info)


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
            data['states'] = StateUpdateFormSet(self.request.POST, instance=self.object)
            data['states'].full_clean()
        else:
            data['states'] = StateUpdateFormSet(instance=self.object)
            if self.object.primary_reference:
                data['form'].fields['reference_doi'].initial = self.object.primary_reference.doi
        return data

    def form_valid(self, form):
        form.instance.updated_by = self.request.user  # login_required in url.py
        return super().form_valid(form)


def protein_table(request):
    ''' renders html for protein table page  '''
    #  return render(request, 'table.html', {"states": State.objects.notdark().select_related('protein').prefetch_related('bleach_measurement')})
    return render(request, 'table.html', {"proteins": Protein.objects.all().prefetch_related('states', 'states__bleach_measurement')})


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


def add_reference(request):
    try:
        slug = request.POST.get('protein')
        doi = request.POST.get('reference_doi')
        P = Protein.objects.get(slug=slug)
        ref, created = Reference.objects.get_or_create(doi=doi)
        P.references.add(ref)
        #P.save()
        return JsonResponse({'url': P.get_absolute_url()})
    except Exception:
        pass


