import datetime

from django.conf import settings
from django.core.cache import cache
from googleapiclient.discovery import build
from oauth2client.service_account import ServiceAccountCredentials


def get_service(
    api_name="analytics",
    api_version="v3",
    scopes=["https://www.googleapis.com/auth/analytics.readonly"],
):
    """Get a service that communicates to a Google API.
    Returns:
        A service that is connected to the specified API.
    """

    keyfile_dict = {
        "type": "service_account",
        "project_id": "fp-base",
        "private_key_id": "0bc1b9014c85ef3799680c4aa1095b57437c2eb4",
        "private_key": settings.GOOGLE_API_PRIVATE_KEY,
        "client_email": settings.GOOGLE_API_CLIENT_EMAIL,
        "client_id": "",
        "auth_uri": "https://accounts.google.com/o/oauth2/auth",
        "token_uri": "https://oauth2.googleapis.com/token",
        "auth_provider_x509_cert_url": "https://www.googleapis.com/oauth2/v1/certs",
        "client_x509_cert_url": (
            "https://www.googleapis.com/robot/v1/metadata/x509/"
            + "626427424458-compute%40developer.gserviceaccount.com"
        ),
    }
    cred = ServiceAccountCredentials.from_json_keyfile_dict(keyfile_dict, scopes=scopes)

    # Build the service object.
    return build(api_name, api_version, credentials=cred)


def get_first_profile_id(service):
    # Use the Analytics service object to get the first profile id.

    # Get a list of all Google Analytics accounts for this user
    accounts = service.management().accounts().list().execute()

    if accounts.get("items"):
        # Get the first Google Analytics account.
        account = accounts.get("items")[0].get("id")

        # Get a list of all the properties for the first account.
        properties = service.management().webproperties().list(accountId=account).execute()

        if properties.get("items"):
            # Get the first property id.
            property = properties.get("items")[0].get("id")

            # Get a list of all views (profiles) for the first property.
            profiles = service.management().profiles().list(accountId=account, webPropertyId=property).execute()

            if profiles.get("items"):
                # return the first view (profile) id.
                return profiles.get("items")[0].get("id")

    return None


def cached_ga_popular(max_age=60 * 60 * 24):
    results = cache.get("ga_popular_proteins")
    if not results:
        service = get_service()

        def f(x):
            return ga_popular_proteins(service, "168069800", x)

        results = {"day": f(1), "week": f(7), "month": f(30), "year": f(365)}
        cache.set("ga_popular_proteins", results, max_age)
    return results


def ga_popular_proteins(service=None, profile_id=None, days=30, max_results=None):
    if service is None:
        service = get_service()
    if not profile_id:
        profile_id = get_first_profile_id(service)
    end_date = datetime.date.today()
    start_date = end_date - datetime.timedelta(days=days)
    data_query = (
        service.data()
        .ga()
        .get(
            **{
                "ids": "ga:" + profile_id,
                "dimensions": "ga:pagePath,ga:pageTitle",
                "metrics": "ga:uniquePageviews",
                "start_date": start_date.strftime("%Y-%m-%d"),
                "end_date": end_date.strftime("%Y-%m-%d"),
                "filters": "ga:pagePath=@protein",
                "sort": "-ga:uniquePageviews",
                "max_results": max_results,
            }
        )
    )
    analytics_data = data_query.execute()
    f = [
        (r[0].replace("/protein/", "").split("/")[0], r[1], r[2])
        for r in analytics_data["rows"]
        if r[0].startswith("/protein")
        and not any([x in r[0] for x in ("bleach", "update")])
        and "not found" not in r[1]
        and " :: " in r[1]
    ]
    total = sum([int(x[2]) for x in f])
    out = [(a, b.split(" ::")[0], 100 * int(c) / total) for a, b, c in f]
    # format percentage: "{:.2%}".format()
    return out
