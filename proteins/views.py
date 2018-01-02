from django.views.generic import DetailView, ListView, CreateView, UpdateView
from django.shortcuts import render, redirect
from django.http import HttpResponseRedirect
from .models import Protein
from references.models import Reference  # breaks application modularity # FIXME
from .forms import ProteinSubmitForm, ProteinUpdateForm, StateFormSet, StateUpdateFormSet
from django.core import serializers
from django.db import transaction


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


class ProteinCreateUpdateMixin:

    def form_valid(self, form):
        # This method is called when valid form data has been POSTed.
        context = self.get_context_data()
        states = context['states']

        with transaction.atomic():
            # only save the form if all the states are also valid
            if states.is_valid():
                self.object = form.save(commit=False)
                doi = form.cleaned_data.get('reference_doi')
                ref, created = Reference.objects.get_or_create(doi=doi)
                self.object.primary_reference = ref

                states.instance = self.object
                saved_states = states.save(commit=False)
                for s in saved_states:
                    if not s.created_by:
                        s.created_by = self.request.user
                    s.updated_by = self.request.user
                    s.save()
                for s in states.deleted_objects:
                    s.delete()

                # FIXME: should allow control of default states in form
                # if only 1 state, make it the default state
                if self.object.states.count() == 1:
                    self.object.default_state = self.object.states.first()
                # otherwise use first non-dark state
                elif self.object.states.count() > 1:
                    for S in self.object.states.all():
                        if not S.is_dark and S not in states.deleted_objects:
                            self.object.default_state = S
                            break
                self.object.save()
            else:
                context.update({
                    'states': states
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
            data['states'] = StateUpdateFormSet(self.request.POST, instance=self.object)
            data['states'].full_clean()
        else:
            data['states'] = StateUpdateFormSet(instance=self.object)
            data['form'].fields['reference_doi'].initial = self.object.primary_reference.doi
        return data

    def form_valid(self, form):
        form.instance.updated_by = self.request.user  # login_required in url.py
        return super().form_valid(form)


def protein_table(request):
    ''' renders html for protein table page  '''
    return render(request, 'table.html', {"proteins": Protein.objects.select_related('default_state')})


def protein_search(request):
    ''' renders html for protein search page  '''
    from .filters import ProteinFilter
    if request.GET:
        f = ProteinFilter(request.GET, queryset=Protein.objects.all().order_by('default_state__em_max').prefetch_related('default_state'))
        print(f.qs)
        if len(f.qs) == 1:
            return redirect(f.qs.first())
    else:
        f = ProteinFilter(request.GET, queryset=Protein.objects.none())
    return render(request, 'proteins/protein_filter.html', {'filter': f})


from django.http import JsonResponse
def add_reference(request):
    try:
        slug = request.POST.get('protein')
        doi = request.POST.get('reference_doi')
        P = Protein.objects.get(slug=slug)
        ref, created = Reference.objects.get_or_create(doi=doi)
        P.references.add(ref)
        P.save()
        return JsonResponse({'url': P.get_absolute_url()})
    except Exception:
        pass

