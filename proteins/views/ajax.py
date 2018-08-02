from django.http import HttpResponseNotAllowed
from django.utils.text import slugify
from django.contrib.auth.decorators import login_required
from django.contrib.admin.views.decorators import staff_member_required
from django.http import JsonResponse
from django.core.mail import mail_managers

from ..models import Protein, Organism, Spectrum, Fluorophore
import reversion


def serialize_comparison(request):
    info = []
    slugs = request.session.get('comparison', [])
    for prot in Protein.objects.filter(slug__in=slugs).prefetch_related('default_state'):
        d = {
            'name': prot.name,
            'slug': prot.slug,
        }
        if prot.default_state:
            d.update({
                'color': prot.default_state.emhex or '#FFFFFF',
                'emMax': prot.default_state.em_max or '',
                'exMax': prot.default_state.ex_max or '',
                'ec': prot.default_state.ext_coeff or '',
                'qy': prot.default_state.qy or '',
                'spectra': prot.d3_spectra(),
            })
        info.append(d)
    return info


def add_to_comparison(request, slug):
    if not request.is_ajax():
        return HttpResponseNotAllowed([])
    if 'comparison' in request.session:
        current = set(request.session['comparison'])
        current.add(slug)
        request.session['comparison'] = list(current)
    else:
        request.session['comparison'] = [slug]
    return JsonResponse({
        'status': 200,
        'comparison_set': serialize_comparison(request)
    })


def remove_from_comparison(request, slug):
    if not request.is_ajax():
        return HttpResponseNotAllowed([])
    if 'comparison' in request.session:
        current = set(request.session['comparison'])
        try:
            current.remove(slug)
        except KeyError:
            pass
        request.session['comparison'] = list(current)
    return JsonResponse({
        'status': 200,
        'comparison_set': serialize_comparison(request)
    })


def get_comparison(request):
    if not request.is_ajax():
        return HttpResponseNotAllowed([])
    return JsonResponse({
        'status': 200,
        'comparison_set': serialize_comparison(request),
    })


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
    if not request.is_ajax():
        return HttpResponseNotAllowed([])

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


def similar_spectrum_owners(request):
    if not request.is_ajax():
        return HttpResponseNotAllowed([])

    name = request.POST.get('owner', None)
    similars = Spectrum.objects.find_similar_owners(name, 0.3)
    data = {
        'similars': [{
            'slug': s.slug,
            'name': s.protein.name if hasattr(s, 'protein') else s.name,
            'url': s.get_absolute_url(),
            'spectra': ([sp.get_subtype_display() for sp in s.spectra.all()]
                         if isinstance(s, Fluorophore)
                         else [s.spectrum.get_subtype_display()])
            }
            for s in similars[:4]]
    }

    return JsonResponse(data)


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
