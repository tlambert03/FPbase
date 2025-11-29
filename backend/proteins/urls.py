from os import getenv

from django.urls import path, re_path
from django.views.generic import TemplateView

from fpbase.decorators import login_required_message_and_redirect as login_required
from proteins import views

app_name = "proteins"

disabled = TemplateView.as_view(template_name="updates_disabled.html")

CONTRIBS_OPEN = not getenv("BLOCK_CONTRIBUTIONS")

urlpatterns = [
    # detail view: /:slug
    re_path(r"^search/", views.protein_search, name="search"),
    re_path(r"^blast/", views.blast_view, name="blast"),
    re_path(
        r"^submit/",
        (
            login_required(
                views.ProteinCreateView.as_view(),
                message="You must be logged in to submit a new protein",
            )
            if CONTRIBS_OPEN
            else disabled
        ),
        name="submit",
    ),
    re_path(
        r"^protein/(?P<slug>[-\w]+)/update/",
        (
            login_required(
                views.ProteinUpdateView.as_view(),
                message="You must be logged in to update protein information",
            )
            if CONTRIBS_OPEN
            else disabled
        ),
        name="update",
    ),
    re_path(r"^table/", views.protein_table, name="table"),
    path(
        "problems/",
        TemplateView.as_view(template_name="problems.html"),
        name="problems",
    ),
    path("problems/gaps/", views.problems_gaps, name="problems-gaps"),
    path(
        "problems/inconsistencies/",
        views.problems_inconsistencies,
        name="problems-inconsistencies",
    ),
    # V2 spectrum submission form (enhanced with client-side processing)
    path(
        "spectra/submit/",
        (
            login_required(
                views.SpectrumCreateViewV2.as_view(),
                message="You must be logged in to submit spectra",
            )
            if CONTRIBS_OPEN
            else disabled
        ),
        name="submit-spectra",
    ),
    # Slug-based submission uses legacy form (pre-selects the protein)
    re_path(
        r"^spectra/submit/(?P<slug>[-\w]+)/$",
        (
            login_required(
                views.SpectrumCreateView.as_view(),
                message="You must be logged in to submit spectra",
            )
            if CONTRIBS_OPEN
            else disabled
        ),
        name="submit-spectra",
    ),
    path(
        "spectra/submitted/",
        views.spectrum_submitted_v2,
        name="spectrum_submitted",
    ),
    # Legacy spectrum submission form
    # NOTE: Must come BEFORE the slug pattern below to avoid "legacy" being matched as a slug
    path(
        "spectra/submit/legacy/",
        (
            login_required(
                views.SpectrumCreateView.as_view(),
                message="You must be logged in to submit a new spectrum",
            )
            if CONTRIBS_OPEN
            else disabled
        ),
        name="submit-spectra-legacy",
    ),
    path(
        "spectra/submitted/legacy/",
        views.spectrum_submitted,
        name="spectrum_submitted_legacy",
    ),
    re_path(
        r"^spectra/submit/legacy/(?P<slug>[-\w]+)/$",
        (
            login_required(
                views.SpectrumCreateView.as_view(),
                message="You must be logged in to submit a new spectrum",
            )
            if CONTRIBS_OPEN
            else disabled
        ),
        name="submit-spectra-legacy",
    ),
    re_path(r"^spectra/(?P<slug>[-\w]+)", views.protein_spectra, name="spectra"),
    path("spectra/", views.protein_spectra, name="spectra"),
    path("spectra-graph/", views.protein_spectra_graph, name="spectra_graph"),
    re_path(r"^spectra_csv/", views.spectra_csv, name="spectra_csv"),
    re_path(
        r"^spectra_img/(?P<slug>[-\w]+)(\.(?P<extension>(png)|(svg)|(tif?f)|(pdf)|(jpe?g))?)?$",
        views.spectra_image,
        name="spectra-img",
    ),
    path(
        "spectra_url_builder/",
        TemplateView.as_view(template_name="spectrum_url_form.html"),
        name="spectra-url-builder",
    ),
    re_path(r"^fret/", views.fret_chart, name="fret"),
    path("compare/", views.ComparisonView.as_view(), name="compare"),
    re_path(
        r"^compare/(?P<proteins>[\w,\-]+)/$",
        views.ComparisonView.as_view(),
        name="compare",
    ),
    re_path(
        r"^lineage/",
        TemplateView.as_view(template_name="lineage.html"),
        name="lineage-list",
    ),
    re_path(r"^chart/", TemplateView.as_view(template_name="ichart.html"), name="ichart"),
    re_path(
        r"^collections/(?P<owner>[\w.@+-]+)/?$",
        views.CollectionList.as_view(),
        name="collections",
    ),
    re_path(r"^collections/", views.CollectionList.as_view(), name="collections"),
    path(
        "collection/<int:pk>/",
        views.CollectionDetail.as_view(),
        name="collection-detail",
    ),
    re_path(
        r"^collection/(?P<pk>\d+)/update/",
        login_required(
            views.CollectionUpdateView.as_view(),
            message="You must be logged in to update collections",
        ),
        name="updatecollection",
    ),
    re_path(
        r"^collection/(?P<pk>\d+)/delete/",
        login_required(
            views.CollectionDeleteView.as_view(),
            message="You must be logged in to delete collections",
        ),
        name="deletecollection",
    ),
    re_path(
        r"^collection/create/",
        login_required(
            views.CollectionCreateView.as_view(),
            message="You must be logged in to create a new collection",
        ),
        name="newcollection",
    ),
    path("organisms/", views.OrganismListView.as_view(), name="organism-list"),
    path(
        "organism/<int:pk>/",
        views.OrganismDetailView.as_view(),
        name="organism-detail",
    ),
    path("activity", views.ActivityView.as_view(), name="activity"),
    re_path(
        r"^protein/(?P<slug>[-\w]+)/$",
        views.ProteinDetailView.as_view(),
        name="protein-detail",
    ),
    re_path(
        r"^protein/(?P<slug>[-\w]+)/bleach/$",
        views.protein_bleach_formsets,
        name="protein-bleach-form",
    ),
    re_path(
        r"^protein/(?P<slug>[-\w]+)/history/$",
        views.protein_history,
        name="protein-history",
    ),
    re_path(
        r"^bleach_comparison/(?P<pk>[-\w]+)/$",
        views.bleach_comparison,
        name="bleach-comparison",
    ),
    path("bleach_comparison/", views.bleach_comparison, name="bleach-comparison"),
    re_path(
        r"^protein/(?P<slug>[-\w]+)/rev/(?P<rev>\d+)$",
        views.ProteinDetailView.as_view(),
        name="protein-detail",
    ),
    re_path(
        r"^protein/(?P<slug>[-\w]+)/ver/(?P<ver>\d+)$",
        views.ProteinDetailView.as_view(),
        name="protein-detail",
    ),
    path(
        "autocomplete-protein/",
        views.ProteinAutocomplete.as_view(),
        name="protein-autocomplete",
    ),
    path(
        "autocomplete-lineage/",
        views.LineageAutocomplete.as_view(),
        name="lineage-autocomplete",
    ),
    path(
        "autocomplete-state/",
        views.StateAutocomplete.as_view(),
        name="state-autocomplete",
    ),
    path(
        "autocomplete-filter/",
        views.FilterAutocomplete.as_view(),
        name="filter-autocomplete",
    ),
    re_path(
        r"^microscope/create/",
        login_required(
            views.MicroscopeCreateView.as_view(),
            message="You must be logged in to create a microscope configuration",
        ),
        name="newmicroscope",
    ),
    re_path(
        r"^microscope/(?P<pk>[-\w]+)/delete/",
        login_required(
            views.MicroscopeDeleteView.as_view(),
            message="You must be logged in to delete microscopes",
        ),
        name="deletemicroscope",
    ),
    re_path(
        r"^embedscope/(?P<pk>[-\w]+)/$",
        views.MicroscopeEmbedView.as_view(),
        name="microscope-embed",
    ),
    re_path(
        r"^microscope/(?P<pk>[-\w]+)/report/$",
        views.ScopeReportView.as_view(),
        name="microscope-report",
    ),
    re_path(
        r"^microscope/(?P<pk>[-\w]+)/$",
        views.MicroscopeDetailView.as_view(),
        name="microscope-detail",
    ),
    re_path(
        r"^microscope/(?P<pk>[-\w]+)/update/",
        login_required(
            views.MicroscopeUpdateView.as_view(),
            message="You must be logged in to update microscopes",
        ),
        name="updatemicroscope",
    ),
    re_path(
        r"^microscopes/(?P<owner>[\w.@+-]+)/?$",
        views.MicroscopeList.as_view(),
        name="microscopes",
    ),
    re_path(r"^microscopes/", views.MicroscopeList.as_view(), name="microscopes"),
    # AJAX
    path("ajax/add_taxonomy/", views.add_organism, name="add_taxonomy"),
    re_path(
        r"^ajax/filter_import/(?P<brand>[-\w]+)$",
        views.filter_import,
        name="filter_import",
    ),
    re_path(
        r"^ajax/add_protein_reference/(?P<slug>[-\w]+)/$",
        views.add_reference,
        name="add_protein_reference",
    ),
    re_path(
        r"^ajax/add_protein_excerpt/(?P<slug>[-\w]+)/$",
        views.add_protein_excerpt,
        name="add_protein_excerpt",
    ),
    re_path(
        r"^ajax/admin_approve_protein/(?P<slug>[-\w]+)/$",
        views.approve_protein,
        name="admin_approve_protein",
    ),
    path(
        "ajax/admin_revert_version/<int:ver>",
        views.revert_version,
        name="admin_revert_version",
    ),
    path(
        "ajax/admin_revert_revision/<int:rev>",
        views.revert_revision,
        name="admin_revert_revision",
    ),
    re_path(
        r"^ajax/update_transitions/(?P<slug>[-\w]+)/$",
        views.update_transitions,
        name="update_transitions",
    ),
    path(
        "ajax/validate_proteinname/",
        views.validate_proteinname,
        name="validate_proteinname",
    ),
    path(
        "ajax/validate_spectrumownername/",
        views.similar_spectrum_owners,
        name="validate_spectrumownername",
    ),
    path(
        "ajax/spectrum_preview/",
        views.spectrum_preview,
        name="spectrum_preview",
    ),
    path(
        "pending-spectra/",
        views.pending_spectra_dashboard,
        name="pending_spectra_dashboard",
    ),
    path(
        "ajax/pending_spectrum_action/",
        views.pending_spectrum_action,
        name="pending_spectrum_action",
    ),
    path(
        "ajax/remove_from_collection/",
        views.collection_remove,
        name="collection-remove",
    ),
    path("ajax/flag_object/", views.flag_object, name="flag_object"),
    path("ajax/add_to_collection/", views.add_to_collection, name="add_to_collection"),
    path("ajax/comparison/", views.update_comparison, name="update-comparison"),
    re_path(r"^ajax/lineage/(?P<slug>[-\w]+)/$", views.get_lineage, name="get-lineage"),
    path("ajax/lineage/", views.get_lineage, name="get-lineage"),
    re_path(
        r"^ajax/org_lineage/(?P<org>[-\w]+)/$",
        views.get_lineage,
        name="get-org-lineage",
    ),
    re_path(
        r"^microscope/(?P<pk>[-\w]+)/report/json/$",
        views.scope_report_json,
        name="scope_report_json",
    ),
    re_path(r"^widget/(?P<slug>[-\w]+)/$", views.Widget.as_view(), name="widget-detail"),
]
