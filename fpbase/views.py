from fpbase.forms import ContactForm
from django.views.generic.edit import FormView
from django.shortcuts import render
from django.conf import settings
from proteins.models import Protein
from sentry_sdk import last_event_id
from django.views.decorators.cache import cache_page


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


def server_error(request, *args, **argv):
    return render(request, "500.html", {
        'sentry_event_id': last_event_id(),
        'sentry_js_dsn': getattr(settings, 'SENTRY_JS_DSN', '')
    }, status=500)


@cache_page(10)
def testview(request):
    import logging
    logger = logging.getLogger(__name__)
    p = Protein.objects.get(name='mNeonGreen')
    logger.info(p)
    return render(request, 'pages/test.html', {'protein': p})


