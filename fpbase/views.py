from fpbase.forms import ContactForm
from django.views.generic.edit import FormView
from django.views.generic import TemplateView
from django.shortcuts import render
from django.conf import settings
from proteins.models import Protein, Spectrum
from sentry_sdk import last_event_id
from django.views.decorators.cache import cache_page
from django.contrib.admin.views.decorators import staff_member_required


class HomeView(TemplateView):
    template_name = "pages/home.html"

    def get_context_data(self):
        data = super().get_context_data()
        data["stats"] = {
            "proteins": Protein.objects.count(),
            "protspectra": Spectrum.objects.exclude(owner_state=None).count(),
        }
        return data


class ContactView(FormView):
    template_name = "pages/contact.html"
    form_class = ContactForm
    success_url = "/thanks/"

    def form_valid(self, form):
        # This method is called when valid form data has been POSTed.
        # It should return an HttpResponse.
        if self.request.recaptcha_is_valid:
            form.send_email()
            return super().form_valid(form)
        return render(self.request, self.template_name, self.get_context_data())


@staff_member_required
def test500(request):
    # Return an "Internal Server Error" 500 response code.
    raise Exception("Make response code 500!")


def server_error(request, *args, **argv):
    return render(
        request,
        "500.html",
        {
            "sentry_event_id": last_event_id(),
            "sentry_dsn": getattr(settings, "SENTRY_DSN", ""),
        },
        status=500,
    )


@cache_page(10)
def testview(request):
    import logging

    logger = logging.getLogger(__name__)
    p = Protein.objects.get(name="mNeonGreen")
    logger.info(p)
    return render(request, "pages/test.html", {"protein": p})
