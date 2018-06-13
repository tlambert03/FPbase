from fpbase.forms import ContactForm
from django.views.generic.edit import FormView
from django.shortcuts import render


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
