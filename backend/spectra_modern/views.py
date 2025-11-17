"""Views for modern spectrum submission without JS-specific dependencies."""

from __future__ import annotations

import contextlib
import io
import logging
from textwrap import dedent

from django.conf import settings
from django.core.mail import EmailMessage
from django.http import JsonResponse
from django.urls import reverse, reverse_lazy
from django.views.generic import CreateView

from fpbase.util import is_ajax, uncache_protein_page
from proteins.models import Protein, Spectrum, State

from .forms import ModernSpectrumForm

logger = logging.getLogger(__name__)


class ModernSpectrumCreateView(CreateView):
    """
    Modern spectrum creation view with clean separation of concerns.

    This view doesn't require any specific JavaScript libraries. All progressive
    enhancement is handled on the client side using data attributes.
    """

    model = Spectrum
    form_class = ModernSpectrumForm
    template_name = "spectra_modern/spectrum_form.html"

    def get_success_url(self, **kwargs):
        return reverse_lazy("spectra_modern:submitted")

    def get_initial(self):
        """Pre-populate form if coming from protein page."""
        init = super().get_initial()
        if self.kwargs.get("slug", False):
            with contextlib.suppress(Exception):
                self.protein = Protein.objects.get(slug=self.kwargs.get("slug"))
                init["owner_state"] = self.protein.default_state
                init["category"] = Spectrum.PROTEIN
        return init

    def get_form_kwargs(self):
        """Pass user to form."""
        kwargs = super().get_form_kwargs()
        kwargs["user"] = self.request.user
        return kwargs

    def get_form(self, form_class=None):
        """Customize form if coming from protein page."""
        form = super().get_form()

        if self.kwargs.get("slug", False):
            with contextlib.suppress(Exception):
                # Limit state choices to this protein
                form.fields["owner_state"].queryset = State.objects.filter(protein=self.protein).select_related(
                    "protein"
                )
                form.fields["owner_state"].required = True
                form.fields["owner_state"].empty_label = None
                form.fields["category"].disabled = True
        return form

    def form_valid(self, form):
        """Handle successful form submission."""
        # Set status to pending for non-staff users
        if not self.request.user.is_staff:
            form.instance.status = Spectrum.STATUS.pending

        response = super().form_valid(form)

        # Clear protein page cache
        with contextlib.suppress(Exception):
            if self.object.owner_state:
                uncache_protein_page(self.object.owner_state.protein.slug, self.request)

        # Send notification email for non-staff submissions
        if not self.request.user.is_staff:
            body = f"""
            A new spectrum has been submitted to FPbase.

            Admin URL: {self.request.build_absolute_uri(form.instance.get_admin_url())}
            User: {self.request.user}
            Data:

            {form.cleaned_data}
            """
            EmailMessage(
                subject=f"[FPbase] Spectrum needs validation: {form.cleaned_data.get('owner', 'Unknown')}",
                body=dedent(body),
                to=[a[1] for a in settings.ADMINS],
                headers={"X-Mailgun-Track": "no"},
            ).send()

        self.request.session["spectrum_name"] = self.object.name
        return response

    def get_context_data(self, **kwargs):
        """Add extra context for template."""
        context = super().get_context_data(**kwargs)
        if hasattr(self, "protein"):
            context["protein"] = self.protein

        # Add URLs for AJAX endpoints as data for progressive enhancement
        context["state_autocomplete_url"] = reverse("proteins:state-autocomplete")
        context["validate_owner_url"] = reverse("proteins:validate_spectrumownername")
        context["preview_url"] = reverse("spectra_modern:preview")

        return context


def spectrum_submitted(request):
    """Simple success page after spectrum submission."""
    from django.shortcuts import render

    context = {"spectrum_name": request.session.get("spectrum_name", "")}
    return render(request, "spectra_modern/spectrum_submitted.html", context)


def spectrum_preview(request) -> JsonResponse:
    """
    AJAX endpoint to preview spectrum data with server-side normalization.

    This is a clean JSON API endpoint that doesn't make assumptions about
    what JavaScript is running on the client.
    """
    if request.method != "POST":
        return JsonResponse({"error": "POST required"}, status=405)

    try:
        logger.debug("Spectrum preview request from user: %s", request.user)

        # Determine data source
        data_source = request.POST.get("data_source", "file")
        use_manual_data = data_source == "manual"
        logger.debug("Data source: %s", data_source)

        # Create form with appropriate data
        if use_manual_data:
            form = ModernSpectrumForm(request.POST, None, user=request.user)
        else:
            post_data = request.POST.copy()
            post_data["data"] = ""
            form = ModernSpectrumForm(post_data, request.FILES, user=request.user)

        if not form.is_valid():
            logger.warning("Form validation failed: %s", form.errors)
            return JsonResponse(
                {
                    "error": "Form validation failed. Please check your input data.",
                    "form_errors": dict(form.errors),
                },
                status=400,
            )

        # Create temporary spectrum instance (don't save to DB)
        temp_spectrum: Spectrum = form.save(commit=False)
        temp_spectrum.created_by = request.user

        file_was_processed = not use_manual_data and bool(request.FILES)

        # Run normalization
        try:
            logger.debug("Running spectrum normalization...")
            temp_spectrum.clean()
            logger.debug("Normalization completed")
        except Exception as e:
            logger.exception("Data processing failed: %s", e)
            return JsonResponse(
                {
                    "error": f"Data processing failed: {e!s}",
                    "details": "The spectrum data could not be normalized.",
                },
                status=400,
            )

        # Generate SVG preview
        try:
            logger.debug("Generating SVG...")
            svg_buffer = io.BytesIO()
            temp_spectrum.spectrum_img(
                fmt="svg",
                output=svg_buffer,
                ylabels=True,
                figsize=(12, 4.5),
            )
            svg_data = svg_buffer.getvalue().decode("utf-8")
            logger.debug("SVG generated successfully")
        except Exception as e:
            logger.exception("SVG generation failed: %s", e)
            return JsonResponse({"error": f"Chart generation failed: {e!s}"}, status=500)

        # Return preview data
        preview_data = {
            "svg": svg_data,
            "peak_wave": temp_spectrum.peak_wave,
            "min_wave": temp_spectrum.min_wave,
            "max_wave": temp_spectrum.max_wave,
            "name": temp_spectrum.name,
            "category": temp_spectrum.category,
            "subtype": temp_spectrum.subtype,
            "data_points": len(temp_spectrum.data),
            "raw_data": temp_spectrum.data,
            "file_was_processed": file_was_processed,
            "normalization_applied": True,
        }

        logger.debug("Preview generated successfully")
        return JsonResponse(
            {
                "success": True,
                "preview": preview_data,
                "message": f"{temp_spectrum.name!r} spectrum processed successfully.",
            }
        )

    except Exception as e:
        logger.exception("Unexpected error in spectrum_preview: %s", e)
        return JsonResponse(
            {
                "error": f"Unexpected server error: {e!s}",
                "traceback": str(e) if request.user.is_staff else None,
            },
            status=500,
        )


def state_search(request) -> JsonResponse:
    """
    Simple JSON endpoint for protein state search.

    This replaces django-autocomplete-light with a clean JSON API
    that any frontend can consume.
    """
    if not is_ajax(request):
        return JsonResponse({"error": "AJAX required"}, status=400)

    query = request.GET.get("q", "").strip()
    if len(query) < 2:
        return JsonResponse({"results": []})

    # Search for states by protein name
    states = (
        State.objects.filter(protein__name__icontains=query)
        .select_related("protein")
        .order_by("protein__name", "name")[:20]
    )

    results = [
        {
            "id": state.id,
            "text": f"{state.protein.name} - {state.name}",
            "protein_name": state.protein.name,
            "state_name": state.name,
        }
        for state in states
    ]

    return JsonResponse({"results": results})
