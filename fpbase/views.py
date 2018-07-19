from fpbase.forms import ContactForm
from django.views.generic.edit import FormView
from django.shortcuts import render
from django.http import HttpResponseServerError
from django.template import loader


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
    return HttpResponseServerError(template.render({'request': request}))
