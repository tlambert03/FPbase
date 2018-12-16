import reversion
from django.views.generic import DetailView, CreateView, UpdateView, ListView, base
from django.views.decorators.cache import cache_page
from django.views.decorators.csrf import csrf_protect
from django.views.decorators.vary import vary_on_cookie
from django.core.mail import mail_managers, mail_admins
from django.shortcuts import render, get_object_or_404
from django.http import (HttpResponseRedirect, HttpResponse, JsonResponse,
                         HttpResponseBadRequest, HttpResponseNotAllowed,
                         Http404)
from django.forms.models import modelformset_factory
from django.contrib.auth.decorators import login_required
from django.contrib.admin.views.decorators import staff_member_required
from django.db import transaction
from django.db.models import Case, When, Count, Prefetch, Q
from django.shortcuts import redirect
from django.apps import apps
from django.utils.html import strip_tags
from django.utils.decorators import method_decorator

from fpbase.util import uncache_protein_page
from ..models import Protein, State, Organism, BleachMeasurement, Spectrum, Excerpt
from ..forms import (ProteinForm, StateFormSet, StateTransitionFormSet, LineageFormSet,
                     BleachMeasurementForm, bleach_items_formset, BleachComparisonForm)

from proteins.util.spectra import spectra2csv
from proteins.util.maintain import check_lineages
from proteins.util.helpers import most_favorited, link_excerpts
from proteins.extrest.ga import cached_ga_popular
from references.models import Reference  # breaks application modularity
from reversion.views import _RollBackRevisionView
from reversion.models import Version
import json
import io
from django.utils.safestring import mark_safe
from django.utils.text import slugify
from django.db.models import Func, F, Value
from django.contrib import messages
import operator
import logging

logger = logging.getLogger(__name__)


class ProteinDetailView(DetailView):
    ''' renders html for single protein page  '''
    queryset = Protein.visible.annotate(has_spectra=Count('states__spectra')) \
        .prefetch_related('states', 'excerpts__reference', 'oser_measurements__reference') \
        .select_related('primary_reference')

    @method_decorator(cache_page(60 * 30))
    @method_decorator(vary_on_cookie)
    def dispatch(self, *args, **kwargs):
        return super().dispatch(*args, **kwargs)

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
        try:
            return super().get(request, *args, **kwargs)
        except Http404:
            name = slugify(self.kwargs.get(self.slug_url_kwarg))
            aliases_lower = Func(Func(F('aliases'), function='unnest'), function='LOWER')
            remove_space = Func(aliases_lower, Value(' '), Value('-'), function='replace')
            final = Func(remove_space, Value('.'), Value(''), function='replace')
            D = dict(Protein.objects.annotate(aka=final).values_list('aka', 'id'))
            if name in D:
                obj = Protein.objects.get(id=D[name])
                messages.add_message(self.request, messages.INFO,
                                     'The URL {} was not found.  You have been forwarded here'
                                     .format(self.request.get_full_path()))
                return HttpResponseRedirect(obj.get_absolute_url())
            raise

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
        data['spectra_ids'] = ",".join([str(sp.id) for state in self.object.states.all() for sp in state.spectra.all()])
        # data['additional_references'] = Reference.objects.filter(proteins=self.object).exclude(id=self.object.primary_reference.id).order_by('-year')

        # put links in excerpts
        data['excerpts'] = link_excerpts(self.object.excerpts.all(), self.object.name, self.object.aliases)
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
        lineage = context['lineage']

        if states.is_valid() and lineage.is_valid():
            # check if another protein already has the sequence that would be
            if not form.cleaned_data.get('seq'):
                for lform in lineage.forms:
                    mut = lform.cleaned_data.get('mutation')
                    par = lform.cleaned_data.get('parent')
                    if par and mut:
                        seq = par.protein.seq.mutate(mut)
                        try:
                            prot = Protein.objects.get(seq__iexact=str(seq))
                            msg = mark_safe('<a href="{}" style="text-decoration: underline;">{}</a> already has the sequence that would be generated by this parent & mutation'.format(
                                            prot.get_absolute_url(), prot.name))
                            lform.add_error(None, msg)
                            context.update({
                                'states': states,
                                'lineage': lineage,
                            })
                            return self.render_to_response(context)
                        except Protein.DoesNotExist:
                            pass

            with transaction.atomic():
                with reversion.create_revision():
                    self.object = form.save()
                    doi = form.cleaned_data.get('reference_doi')
                    if doi:
                        ref, created = Reference.objects.get_or_create(doi=doi.lower())
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

                    lineage.instance = self.object
                    for lin in lineage.save(commit=False):
                        # if the form has been cleared and there are no children,
                        # let's clean up a little
                        if not (lin.mutation and lin.parent) and not lin.children.exists():
                            lin.delete()
                        else:
                            if not lin.created_by:
                                lin.created_by = self.request.user
                            lin.updated_by = self.request.user
                            lin.reference = self.object.primary_reference
                            lin.save()
                    for lin in lineage.deleted_objects:
                        lin.delete()

                    if hasattr(self.object, 'lineage'):
                        if not self.object.seq:
                            seq = self.object.lineage.parent.protein.seq.mutate(self.object.lineage.mutation)
                            self.object.seq = str(seq)
                        if not self.object.parent_organism:
                            self.object.parent_organism = self.object.lineage.root_node.protein.parent_organism

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
                    try:
                        uncache_protein_page(self.object.slug, self.request)
                    except Exception as e:
                        logger.error('failed to uncache protein: {}'.format(e))
        else:
            context.update({
                'states': states,
                'lineage': lineage,
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
            data['lineage'] = LineageFormSet(self.request.POST)
        else:
            data['states'] = StateFormSet()
            data['lineage'] = LineageFormSet()
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

            data['lineage'] = LineageFormSet(
                self.request.POST, instance=self.object)
        else:
            data['states'] = StateFormSet(instance=self.object)
            data['lineage'] = LineageFormSet(instance=self.object)
            if self.object.primary_reference:
                data['form'].fields['reference_doi'].initial = self.object.primary_reference.doi
        return data

    def form_valid(self, form):
        form.instance.updated_by = self.request.user
        return super().form_valid(form)


class ActivityView(ListView):
    template_name = "proteins/activity.html"
    stateprefetch = Prefetch('states', queryset=State.objects.prefetch_related('spectra').order_by('-is_dark', 'em_max'))
    queryset = Protein.visible.prefetch_related(stateprefetch, 'primary_reference').order_by('-created')[:18]

    def get_context_data(self, **kwargs):
        data = super().get_context_data(**kwargs)
        stateprefetch = Prefetch('states', queryset=State.objects.prefetch_related('spectra').order_by('-is_dark', 'em_max'))
        data['proteins_by_date'] = Protein.visible.annotate(nstates=Count('states')) \
            .prefetch_related(stateprefetch, 'primary_reference') \
            .order_by(F('primary_reference__date').desc(nulls_last=True))[:15]
        data['most_popular'] = cached_ga_popular()
        data['most_favorited'] = most_favorited(max_results=12)
        return data


@cache_page(60 * 10)
def spectra_image(request, slug, **kwargs):
    protein = get_object_or_404(Protein.objects.select_related('default_state'), slug=slug)
    try:
        D = {}
        for k, v in request.GET.dict().items():
            if k == 'xlim':
                tmp = [int(x) for x in v.split(',')][:2]
                if len(tmp) == 1:
                    tmp.append(750)
                D[k] = tmp
            else:
                if v.lower() in ('false', 'no'):
                    D[k] = False
                elif v.lower() in ('true', 'yes'):
                    D[k] = True
                else:
                    D[k] = int(v) if v.isdigit() else v
    except Exception as e:
        return HttpResponseBadRequest('failed to parse url parameters as JSON')
    try:
        fmt = kwargs.get('extension', 'png')
        byt = protein.spectra_img(fmt, output=io.BytesIO(), **D)
    except Exception as e:
        logger.error(e)
        return HttpResponseBadRequest(
            'failed to parse url parameters as JSON: {}'.format(e))
    if byt:
        byt.seek(0)
        if fmt == 'svg':
            fmt += '+xml'
        return HttpResponse(byt, content_type="image/%s" % fmt)
    raise Http404()


@cache_page(60 * 10)
@csrf_protect
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


class ComparisonView(base.TemplateView):

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
        prots_with_seqs = prots.exclude(seq__isnull=True)
        if prots_with_seqs.count() > 2:
            a, _ = prots_with_seqs.to_tree('html')
            context['alignment'] = "\n".join(a.splitlines()[3:]).replace("FFEEE0", "FFFFFF")
        elif prots_with_seqs.count() == 2:
            seqA = prots_with_seqs[0]
            seqB = prots_with_seqs[1]
            if seqA.seq and seqB.seq:
                algn = seqA.seq.align_to(seqB.seq)
                out = []
                for i, row in enumerate(str(algn).splitlines()):
                    head = ''
                    if i % 4 == 0:
                        head = prots[0].name
                    elif i % 2 == 0:
                        head = prots[1].name
                    out.append("{:<18.16}{}".format(head if len(head) < 17 else head[:13] + '...', row))
                out.append('\n')
                context['alignment'] = "\n".join(out)
                context['mutations'] = str(seqA.seq.mutations_to(seqB.seq))
        else:
            context['alignment'] = None
        prots = list(prots)
        context['proteins'] = prots
        refs = Reference.objects.filter(primary_proteins__in=prots)
        refs = refs | Reference.objects.filter(proteins__in=prots)
        refs = refs.distinct('id').order_by('id')
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
    linerrors = check_lineages()[0]

    from functools import reduce
    query = reduce(operator.or_, (Q(primary_reference__title__icontains = item)
                   for item in ['activat', 'switch', 'convert', 'dark', 'revers']))

    return render(
        request,
        'seq_problems.html',
        {
            "histags": Protein.objects.filter(seq__icontains='HHHHH').values('name', 'slug'),
            "noseqs": Protein.objects.filter(seq__isnull=True).values('name', 'slug'),
            "nostates": Protein.objects.filter(states=None).values('name', 'slug'),
            "linprobs": [(node.protein, v) for node, v in linerrors.items()],
            "nomet": Protein.objects.exclude(seq__isnull=True).exclude(seq__istartswith='M'),
            "noparent": Protein.objects.filter(parent_organism__isnull=True),
            "only2p": (State.objects.filter(spectra__subtype='2p')
                                    .exclude(spectra__subtype='ex')
                                    .distinct('protein')
                                    .values('protein__name', 'protein__slug')),
            'nolineage': Protein.objects.filter(lineage=None).annotate(ns=Count('states__spectra')).order_by('-ns'),
            'switchers': Protein.objects.annotate(ns=Count('states')).filter(query).filter(ns=1),
            'request': request
        })


def add_reference(request, slug=None):
    try:
        with reversion.create_revision():
            doi = request.POST.get('reference_doi').lower()
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
            try:
                uncache_protein_page(slug, request)
            except Exception as e:
                logger.error('failed to uncache protein: {}'.format(e))
        return JsonResponse({'status': 'success'})
    except Exception as e:
        return JsonResponse({'status': 'failed', 'msg': e})


def add_excerpt(request, slug=None):
    try:
        with reversion.create_revision():
            doi = request.POST.get('excerpt_doi').lower()
            P = Protein.objects.get(slug=slug)
            content = request.POST.get('excerpt_content')
            if content:
                ref, created = Reference.objects.get_or_create(doi=doi)
                P.references.add(ref)
                Excerpt.objects.create(protein=P, reference=ref, content=strip_tags(content), created_by=request.user)
                if not request.user.is_staff:
                    msg = "User: {}\nProtein: {}\nReference: {}, {}\nExcerpt: {}\n{}".format(
                        request.user.username, P, ref, ref.title, strip_tags(content),
                        request.build_absolute_uri(P.get_absolute_url()))
                    mail_managers('Excerpt Added', msg, fail_silently=True)
                P.save()
                reversion.set_user(request.user)
                reversion.set_comment('Excerpt from {} added to {}'.format(ref, P))
                try:
                    uncache_protein_page(slug, request)
                except Exception as e:
                    logger.error('failed to uncache protein: {}'.format(e))
        return JsonResponse({'status': 'success'})
    except Exception as e:
        return JsonResponse({'status': 'failed', 'msg': e})


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
                    try:
                        uncache_protein_page(slug, request)
                    except Exception as e:
                        logger.error('failed to uncache protein: {}'.format(e))
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
                    try:
                        uncache_protein_page(slug, request)
                    except Exception as e:
                        logger.error('failed to uncache protein: {}'.format(e))
            return HttpResponseRedirect(protein.get_absolute_url())
        else:
            return render(request, template_name, {'formset': formset, 'protein': protein})
    else:
        formset = BleachMeasurementFormSet(queryset=qs)
        formset.form.base_fields['state'].queryset = State.objects.filter(protein__slug=slug)
    return render(request, template_name, {'formset': formset, 'protein': protein})


@login_required
def bleach_comparison(request, pk=None):
    template_name = 'proteins/bleach_comparison_form.html'
    if request.method == 'POST':
        formset = bleach_items_formset(request.POST)
        bcf = BleachComparisonForm(request.POST)
        if formset.is_valid() and bcf.is_valid():
            D = bcf.cleaned_data
            D['reference'], _ = Reference.objects.get_or_create(doi=D.pop('reference_doi'))
            for form in formset.forms:
                if form.has_changed():
                    D['state'] = form.cleaned_data['state']
                    D['rate'] = form.cleaned_data['rate']
                    BleachMeasurement.objects.create(**D)
            return redirect(D['reference'])
        else:
            return render(request, template_name, {'formset': formset, 'mainform': bcf})
    else:
        formset = bleach_items_formset()
        if pk:
            reference = get_object_or_404(Reference, id=pk)
            bcf = BleachComparisonForm({'reference_doi': reference.doi})
            bcf.fields['reference_doi'].widget.attrs['readonly'] = True
        else:
            bcf = BleachComparisonForm()
    return render(request, template_name, {'formset': formset, 'mainform': bcf})


class OrganismListView(ListView):
    ''' renders html for single reference page  '''
    queryset = Organism.objects.annotate(num_prot=Count('proteins'))


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


@login_required
def flag_object(request):
    if not request.is_ajax():
        return HttpResponseNotAllowed([])
    try:
        model_type = request.POST["target_model"]
        id = request.POST["target_id"]
        model = apps.get_model(model_type)
        obj = get_object_or_404(model, id=id)

        status = None
        if request.POST["flagged"] == '1':
            if request.user.is_staff:
                obj.status = model.STATUS.approved
                status = 'unflagged'
        else:
            obj.status = model.STATUS.flagged
            status = 'flagged'
        obj.save()

        if status:
            try:
                uncache_protein_page(obj.protein.slug, request)
            except Exception as e:
                logger.error('failed to uncache protein: {}'.format(e))

            mail_admins('%s %s' % (model_type, status),
                        "User: {}\nObject: {}\nID: {}\n{}".format(
                            request.user.username, obj, obj.pk, request.build_absolute_uri(obj.get_absolute_url())),
                        fail_silently=True)

        return JsonResponse({
            'status': status,
        })
    except (KeyError, ValueError):
        return HttpResponseBadRequest()
