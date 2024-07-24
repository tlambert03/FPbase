from contextlib import suppress

from django.conf import settings
from django.core.cache import cache
from google.analytics.data_v1beta import BetaAnalyticsDataClient
from google.analytics.data_v1beta.types import (
    DateRange,
    Dimension,
    Filter,
    FilterExpression,
    Metric,
    RunReportRequest,
)
from google.oauth2.service_account import Credentials

from proteins.models import Protein

PROPERTY_ID = "255212585"


def get_client() -> "BetaAnalyticsDataClient":
    """Get a service that communicates to a Google API."""

    keyfile_dict = {
        "type": "service_account",
        "project_id": "fp-base",
        "private_key_id": settings.GOOGLE_API_PRIVATE_KEY_ID,
        "private_key": settings.GOOGLE_API_PRIVATE_KEY,
        "client_email": settings.GOOGLE_API_CLIENT_EMAIL,
        "client_id": "",
        "auth_uri": "https://accounts.google.com/o/oauth2/auth",
        "token_uri": "https://oauth2.googleapis.com/token",
        "auth_provider_x509_cert_url": "https://www.googleapis.com/oauth2/v1/certs",
        "client_x509_cert_url": (
            "https://www.googleapis.com/robot/v1/metadata/x509/" + settings.GOOGLE_API_CLIENT_EMAIL.replace("@", "%40")
        ),
        "universe_domain": "googleapis.com",
    }
    scopes = ["https://www.googleapis.com/auth/analytics.readonly"]
    credentials = Credentials.from_service_account_info(keyfile_dict, scopes=scopes)
    return BetaAnalyticsDataClient(credentials=credentials)


def cached_ga_popular(max_age=60 * 60 * 24):
    results = cache.get("ga_popular_proteins")
    if not results:
        client = get_client()
        results = {
            "day": ga_popular_proteins(client, 1),
            "week": ga_popular_proteins(client, 7),
            "month": ga_popular_proteins(client, 30),
            "year": ga_popular_proteins(client, 365),
        }
        cache.set("ga_popular_proteins", results, max_age)
    return results


def ga_popular_proteins(client: BetaAnalyticsDataClient, days: int = 30) -> list[tuple[str, str, float]]:
    """Return a list of proteins with their page views in the last `days` days.

    Returns a list of tuples, each containing: `(protein slug, protein name, view percentage)`
    Sorted by view count in descending order.
    """
    slug2name: dict[str, str] = {}
    uuid2slug: dict[str, str] = {}
    for item in Protein.objects.all().values("slug", "name", "uuid"):
        uuid2slug[item["uuid"]] = item["slug"]
        slug2name[item["slug"]] = item["name"]
    request = RunReportRequest(
        property=f"properties/{PROPERTY_ID}",
        date_ranges=[DateRange(start_date=f"{days}daysAgo", end_date="today")],
        dimensions=[Dimension(name="pagePath"), Dimension(name="pageTitle")],
        metrics=[Metric(name="screenPageViews")],
        dimension_filter=FilterExpression(
            filter=Filter(
                field_name="pagePath",
                string_filter=Filter.StringFilter(
                    match_type=Filter.StringFilter.MatchType.FULL_REGEXP,
                    value=r"^/protein/[^/]+/$",
                    case_sensitive=False,
                ),
            )
        ),
    )
    response = client.run_report(request)

    slug2count: dict[str, int] = {}
    for row in response.rows:
        with suppress(Exception):
            page_title = row.dimension_values[1].value
            if "not found" in page_title.lower():
                continue
            slug = row.dimension_values[0].value.replace("/protein/", "").split("/")[0]
            slug = uuid2slug.get(slug, slug)
            count = slug2count.setdefault(slug, 0)
            slug2count[slug] = count + int(row.metric_values[0].value)

    total_views = sum(slug2count.values())
    with_percent = sorted(
        ((slug, slug2name.get(slug, slug), 100 * count / total_views) for slug, count in slug2count.items()),
        key=lambda x: x[2],
        reverse=True,
    )

    return with_percent
