import contextlib
import io
import json
import logging
import traceback
from textwrap import dedent

# from django.views.decorators.cache import cache_page
# from django.views.decorators.vary import vary_on_cookie
from django import forms
from django.conf import settings
from django.contrib.admin.views.decorators import staff_member_required
from django.core.mail import EmailMessage
from django.http import Http404, HttpResponse, JsonResponse
from django.shortcuts import render
from django.template.defaultfilters import slugify
from django.urls import reverse_lazy
from django.views.decorators.clickjacking import xframe_options_exempt
from django.views.decorators.http import require_POST
from django.views.generic import CreateView

from fpbase.util import is_ajax, uncache_protein_page

from ..forms import SpectrumForm
from ..models import Filter, Protein, Spectrum, State
from ..util.importers import add_filter_to_database
from ..util.spectra import spectra2csv


# @cache_page(60 * 10)
# @vary_on_cookie
def protein_spectra(request, slug=None):
    """renders html for protein spectra page"""
    template = "spectra.html"

    if is_ajax(request) and slug is not None:
        """slug represents the slug of the OWNER of the spectra...
        for instance, a protein state.
        function returns an array containing ALL of the spectra
        belonging to that owner
        """
        owner = Spectrum.objects.filter_owner(slug)
        if owner.count():
            return JsonResponse({"spectra": [s.d3dict() for s in owner]})
        raise Http404

    return render(request, template)


@xframe_options_exempt
def protein_spectra_graph(request, slug=None):
    """Simple variant of protein_spectra only returns the graph.

    Used for embedding in iframes.
    """
    return render(request, "spectra_graph.html")


def spectrum_submitted(request):
    context = {"spectrum_name": request.session.get("spectrum_name", "")}
    return render(request, "spectrum_submitted.html", context)


class SpectrumCreateView(CreateView):
    model = Spectrum
    form_class = SpectrumForm

    def get_success_url(self, **kwargs):
        return reverse_lazy("proteins:spectrum_submitted")

    def get_initial(self):
        init = super().get_initial()
        if self.kwargs.get("slug", False):
            with contextlib.suppress(Exception):
                self.protein = Protein.objects.get(slug=self.kwargs.get("slug"))
                init["owner_state"] = self.protein.default_state
                init["category"] = Spectrum.PROTEIN

        return init

    def get_form_kwargs(self):
        kwargs = super().get_form_kwargs()
        kwargs["user"] = self.request.user
        return kwargs

    def get_form(self, form_class=None):
        form = super().get_form()

        if self.kwargs.get("slug", False):
            with contextlib.suppress(Exception):
                form.fields["owner_state"] = forms.ModelChoiceField(
                    required=True,
                    label="Protein (state)",
                    empty_label=None,
                    queryset=State.objects.filter(protein=self.protein).select_related("protein"),
                )
                form.fields["category"].disabled = True
        return form

    def form_valid(self, form):
        # Set the status to "pending" before saving the form
        if not self.request.user.is_staff:
            form.instance.status = Spectrum.STATUS.pending
        # This method is called when valid form data has been POSTed.
        # It should return an HttpResponse.
        response = super().form_valid(form)
        with contextlib.suppress(Exception):
            uncache_protein_page(self.object.owner_state.protein.slug, self.request)

        if not self.request.user.is_staff:
            body = f"""
            A new spectrum has been submitted to FPbase.

            Admin URL: {self.request.build_absolute_uri(form.instance.get_admin_url())}
            User: {self.request.user}
            Data:

            {form.cleaned_data}
            """
            EmailMessage(
                subject=f"[FPbase] Spectrum needs validation: {form.cleaned_data['owner']}",
                body=dedent(body),
                to=[a[1] for a in settings.ADMINS],
                headers={"X-Mailgun-Track": "no"},
            ).send()
        self.request.session["spectrum_name"] = self.object.name
        return response

    def get_context_data(self, **kwargs):
        # Call the base implementation first to get a context
        context = super().get_context_data(**kwargs)
        if hasattr(self, "protein"):
            context["protein"] = self.protein
        return context


def spectra_csv(request):
    try:
        idlist = [int(x) for x in request.GET.get("q", "").split(",") if x]
        spectralist = Spectrum.objects.filter(id__in=idlist).select_related(
            "owner_state__protein",
            "owner_dye",
            "owner_filter",
            "owner_light",
            "owner_camera",
        )
        if spectralist:
            return spectra2csv(spectralist)
    except Exception:
        return HttpResponse("malformed spectra csv request")
    else:
        return HttpResponse("malformed spectra csv request")


def spectrum_preview(request) -> JsonResponse:
    """
    AJAX endpoint to preview spectrum data with server-side normalization
    before final submission.
    """
    logger = logging.getLogger(__name__)
    if request.method != "POST":
        return JsonResponse({"error": "POST required"}, status=405)

    try:
        # Log the request for debugging
        logger.debug("Spectrum preview request from user: %s", request.user)
        # Determine data source based on explicit tab selection
        data_source = request.POST.get("data_source", "file")  # Default to file tab
        use_manual_data = data_source == "manual"
        logger.debug("Tab selection - data_source: %s", data_source)

        if use_manual_data:
            # Manual data tab selected - ignore any files, only use manual data
            logger.debug("Using manual data entry (files ignored)")
            manual_data = request.POST.get("data", "")
            if manual_data:
                logger.debug("Manual data submitted: %s...", manual_data[:100])
            form = SpectrumForm(request.POST, None, user=request.user)  # No files
        else:
            # File upload tab selected - ignore manual data, only use files
            logger.debug("Using file upload (manual data ignored)")
            # Clear manual data to avoid any confusion
            post_data = request.POST.copy()
            post_data["data"] = ""
            form = SpectrumForm(post_data, request.FILES, user=request.user)

        if not form.is_valid():
            logger.warning("Form validation failed: %s", form.errors)
            return JsonResponse(
                {
                    "error": "Form validation failed. Please check your input data.",
                    "form_errors": dict(form.errors),
                    "details": "See form_errors for specific field issues",
                },
                status=400,
            )

        # Create a temporary spectrum instance (don't save to DB)
        temp_spectrum: Spectrum = form.save(commit=False)
        temp_spectrum.created_by = request.user

        # Track if we processed a file upload (only possible when file tab is active)
        file_was_processed = not use_manual_data and bool(request.FILES)

        # Run the normalization/cleaning process without saving
        try:
            logger.debug("Running spectrum normalization...")
            temp_spectrum.clean()  # This runs the normalization
            logger.debug("Spectrum normalization completed successfully")
        except Exception as e:
            logger.exception("Data processing failed: %s", e)
            return JsonResponse(
                {
                    "error": f"Data processing failed: {e!s}",
                    "details": "The spectrum data could not be normalized. Please check your data format.",
                },
                status=400,
            )

        # Generate SVG image using existing matplotlib renderer
        try:
            logger.debug("Generating SVG image...")
            # Use custom parameters for preview: Y-axis labels, proper sizing for web display
            svg_buffer = io.BytesIO()
            temp_spectrum.spectrum_img(
                fmt="svg",
                output=svg_buffer,
                ylabels=True,
                figsize=(12, 4.5),  # Wider to accommodate Y-axis labels
            )
            svg_data = svg_buffer.getvalue().decode("utf-8")
            logger.debug("SVG generation completed successfully")
        except Exception as e:
            logger.exception("SVG generation failed: %s", e)
            return JsonResponse(
                {
                    "error": f"Chart generation failed: {e!s}",
                    "details": "Could not generate spectrum visualization",
                },
                status=500,
            )

        # Generate preview data for display
        try:
            preview_data = {
                "svg": svg_data,
                "peak_wave": temp_spectrum.peak_wave,
                "min_wave": temp_spectrum.min_wave,
                "max_wave": temp_spectrum.max_wave,
                "name": temp_spectrum.name,
                "category": temp_spectrum.category,
                "subtype": temp_spectrum.subtype,
                "data_points": len(temp_spectrum.data),
                "raw_data": temp_spectrum.data,  # Include raw data for editing
                "file_was_processed": file_was_processed,  # Tell frontend if file was processed
                "normalization_applied": True,  # We always normalize in clean()
            }

            logger.debug("Spectrum preview generated successfully")
            return JsonResponse(
                {
                    "success": True,
                    "preview": preview_data,
                    "message": f"{temp_spectrum.name!r} spectrum processed successfully.",
                }
            )

        except Exception as e:
            logger.exception("Preview data generation failed: %s", e)
            return JsonResponse(
                {
                    "error": f"Preview data generation failed: {e!s}",
                    "details": "Could not extract spectrum properties",
                },
                status=500,
            )

    except Exception as e:
        logger.exception("Unexpected error in spectrum_preview: %s", e)
        return JsonResponse(
            {
                "error": f"Unexpected server error: {e!s}",
                "details": "An unexpected error occurred. Please try again or contact support.",
                "traceback": traceback.format_exc() if request.user.is_staff else None,
            },
            status=500,
        )


def filter_import(request, brand):
    part = request.POST["part"]
    new_objects = []
    errors = []
    response = {"status": 0}

    with contextlib.suppress(Filter.DoesNotExist):
        Filter.objects.get(slug=slugify(f"{brand} {part}"))
        response["message"] = f"{part} is already in the database"
        return JsonResponse(response)

    try:
        new_objects, errors = add_filter_to_database(brand, part, request.user)
    except Exception as e:
        response["message"] = str(e)

    if new_objects:
        spectrum = new_objects[0]
        response = {
            "status": 1,
            "objects": spectrum.name,
            "spectra_options": json.dumps(
                {
                    "category": spectrum.category,
                    "subtype": spectrum.subtype,
                    "slug": spectrum.owner.slug,
                    "name": spectrum.owner.name,
                }
            ),
        }
    elif errors:
        with contextlib.suppress(Exception):
            if errors[0][1].as_data()["owner"][0].code == "owner_exists":
                response["message"] = f"{part} already appears to be imported"
    return JsonResponse(response)


@staff_member_required
def pending_spectra_dashboard(request):
    """Admin dashboard for reviewing pending spectra submissions."""
    pending_spectra = (
        Spectrum.objects.all_objects()
        .filter(status=Spectrum.STATUS.pending)
        .select_related(
            "created_by",
            "owner_state__protein",
            "owner_dye",
            "owner_filter",
            "owner_camera",
            "owner_light",
        )
        .order_by("-created")
    )

    # Generate preview SVGs for each spectrum
    spectra_data = []
    for spectrum in pending_spectra:
        svg_preview = None
        if spectrum.data is not None:
            try:
                buffer = spectrum.spectrum_img(fmt="svg", xlabels=True, ylabels=True, figsize=(10, 3))
                svg_preview = buffer.getvalue().decode("utf-8")
            except Exception:
                svg_preview = None

        spectra_data.append(
            {
                "id": spectrum.id,
                "name": spectrum.name,
                "owner": spectrum.owner,
                "owner_type": spectrum.get_category_display(),
                "category": spectrum.get_category_display(),
                "subtype": spectrum.get_subtype_display(),
                "source": spectrum.source or "N/A",
                "created": spectrum.created,
                "created_by": spectrum.created_by,
                "created_by_email": spectrum.created_by.email if spectrum.created_by else None,
                "svg_preview": svg_preview,
                "admin_url": f"/admin/proteins/spectrum/{spectrum.id}/change/",
            }
        )

    context = {
        "spectra": spectra_data,
        "count": len(spectra_data),
    }

    return render(request, "pending_spectra_dashboard.html", context)


@staff_member_required
@require_POST
def pending_spectrum_action(request):
    """Handle actions (accept/reject/delete) on pending spectra."""
    try:
        spectrum_ids = request.POST.getlist("spectrum_ids[]")
        action = request.POST.get("action")

        if not spectrum_ids or not action:
            return JsonResponse({"success": False, "error": "Missing spectrum_ids or action"}, status=400)

        # For revert (undo), we need to find spectra regardless of status
        if action == "revert":
            spectra = Spectrum.objects.all_objects().filter(id__in=spectrum_ids)
            if not spectra.exists():
                return JsonResponse({"success": False, "error": "No spectra found with provided IDs"}, status=404)
            count = spectra.count()
            spectra.update(status=Spectrum.STATUS.pending)
            message = f"Reverted {count} spectrum(s) to pending"
        else:
            # For other actions, only work with pending spectra
            spectra = Spectrum.objects.all_objects().filter(id__in=spectrum_ids, status=Spectrum.STATUS.pending)

            if not spectra.exists():
                return JsonResponse({"success": False, "error": "No pending spectra found with provided IDs"}, status=404)

            count = spectra.count()

            if action == "accept":
                spectra.update(status=Spectrum.STATUS.approved)
                # Clear cache for affected protein pages
                for spectrum in spectra:
                    with contextlib.suppress(Exception):
                        if spectrum.owner_state:
                            uncache_protein_page(spectrum.owner_state.protein.slug, request)
                message = f"Accepted {count} spectrum(s)"

            elif action == "reject":
                spectra.update(status=Spectrum.STATUS.rejected)
                message = f"Rejected {count} spectrum(s)"

            elif action == "delete":
                spectra.delete()
                message = f"Deleted {count} spectrum(s)"

            else:
                return JsonResponse({"success": False, "error": f"Unknown action: {action}"}, status=400)

        return JsonResponse({"success": True, "message": message, "count": count})

    except Exception as e:
        logger = logging.getLogger(__name__)
        logger.exception("Error in pending_spectrum_action: %s", e)
        return JsonResponse({"success": False, "error": str(e)}, status=500)


@staff_member_required
@require_POST
def pending_spectrum_email(request):
    """Send email to submitter(s) of pending spectra."""
    try:
        spectrum_ids = request.POST.getlist("spectrum_ids[]")
        email_subject = request.POST.get("subject", "")
        email_body = request.POST.get("body", "")

        if not spectrum_ids or not email_subject or not email_body:
            return JsonResponse(
                {"success": False, "error": "Missing spectrum_ids, subject, or body"},
                status=400,
            )

        spectra = Spectrum.objects.all_objects().filter(id__in=spectrum_ids, status=Spectrum.STATUS.pending)

        if not spectra.exists():
            return JsonResponse({"success": False, "error": "No pending spectra found with provided IDs"}, status=404)

        # Get unique submitters
        submitters = {}
        for spectrum in spectra:
            if spectrum.created_by and spectrum.created_by.email:
                if spectrum.created_by.email not in submitters:
                    submitters[spectrum.created_by.email] = {
                        "user": spectrum.created_by,
                        "spectra": [],
                    }
                submitters[spectrum.created_by.email]["spectra"].append(spectrum)

        if not submitters:
            return JsonResponse({"success": False, "error": "No submitters with email addresses found"}, status=400)

        # Send individual emails to each submitter
        emails_sent = 0
        for email_address, data in submitters.items():
            user = data["user"]
            user_spectra = data["spectra"]

            # Personalize email body with spectrum list
            spectrum_list = "\n".join([f"- {s.name} (ID: {s.id})" for s in user_spectra])
            personalized_body = f"""Hi {user.username},

{email_body}

Spectrum(s) in question:
{spectrum_list}

Best regards,
The FPbase Team
"""

            EmailMessage(
                subject=f"[FPbase] {email_subject}",
                body=personalized_body,
                to=[email_address],
                reply_to=[a[1] for a in settings.ADMINS],
                headers={"X-Mailgun-Track": "no"},
            ).send()
            emails_sent += 1

        return JsonResponse(
            {
                "success": True,
                "message": f"Sent {emails_sent} email(s) to {emails_sent} submitter(s)",
                "emails_sent": emails_sent,
            }
        )

    except Exception as e:
        logger = logging.getLogger(__name__)
        logger.exception("Error in pending_spectrum_email: %s", e)
        return JsonResponse({"success": False, "error": str(e)}, status=500)
