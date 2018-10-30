from fpbase.forms import ContactForm
from django.views.generic.edit import FormView
from django.shortcuts import render
from django.http import HttpResponseServerError
from django.template import loader
from django.conf import settings


class ContactView(FormView):
    template_name = 'pages/contact.html'
    form_class = ContactForm
    success_url = '/thanks/'

    def form_valid(self, form):
        # This method is called when valid form data has been POSTed.
        # It should return an HttpResponse.
        if self.request.recaptcha_is_valid:
            form.send_email()
            return super().form_valid(form)
        return render(self.request, self.template_name, self.get_context_data())


def test500(request):
        # Return an "Internal Server Error" 500 response code.
        raise Exception('Make response code 500!')


def server_error(request, template_name='500.html'):
    template = loader.get_template(template_name)
    return HttpResponseServerError(template.render({
        'request': request,
        'sentry_js_dsn': getattr(settings, 'SENTRY_JS_DSN', '')
    }))


import json
from mptt.templatetags.mptt_tags import cache_tree_children
from proteins.models import Protein, Lineage
from django.views.generic import TemplateView
from collections import defaultdict

def recursive_node_to_dict(node, widths=None):
    if not widths:
        widths = defaultdict(int)
    widths[node.level] += 1

    result = {
        'name': node.protein.name,
        'mut': ", ".join(node.mutation.split('/')),
        # 'mut': node.display_mutation(maxwidth=10) or "null",
        # 'url': node.protein.get_absolute_url(),
        'url': '/test/' + node.protein.slug,
        'bg': node.protein.emhex,
    }

    children = []
    for c in node.get_children():
        child, widths = recursive_node_to_dict(c, widths)
        children.append(child)
    if children:
        result['children'] = children
    return result, widths


class testview(TemplateView):
    template_name = 'pages/test.html'

    def get_context_data(self, slug=None):
        data = super().get_context_data()

        if slug:
            item = Lineage.objects.get(protein__slug=slug)
            ids = item.get_family()
        else:
            ids = Lineage.objects.all().values_list('id', flat=True)
        # cache upfront everything we're going to need
        root_nodes = Lineage.objects\
            .filter(id__in=ids)\
            .select_related('protein', 'protein__default_state')\
            .prefetch_related('protein__states')\
            .get_cached_trees()

        D = {
            'name': 'fakeroot',
            'children': [],
            'widths': defaultdict(int)
        }
        for n in root_nodes:
            result, D['widths'] = recursive_node_to_dict(n, D['widths'])
            if 'children' in result:
                D['children'].append(result)
        D['max_width'] = max(D['widths'].values())
        D['max_depth'] = max(D['widths'])
        data['tree'] = json.dumps(D)
        return data


