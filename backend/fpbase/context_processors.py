from django.conf import settings


def api_keys(request):
    return {
        "sentry_dsn": getattr(settings, "SENTRY_DSN", ""),
        "algolia_suffix": getattr(settings, "ALGOLIA_SUFFIX", "prod"),
        "algolia_public_key": getattr(settings, "ALGOLIA_PUBLIC_KEY", ""),
        "algolia_app_id": getattr(settings, "ALGOLIA", {}).get("APPLICATION_ID"),
    }


def canonical(request):
    sha = settings.HEROKU_SLUG_COMMIT
    return {
        "CANONICAL_URL": settings.CANONICAL_URL,
        "HELP_URL": "https://help.fpbase.org/",
        "ABSOLUTE_ROOT": request.build_absolute_uri("/")[:-1].strip("/"),
        "ABSOLUTE_ROOT_URL": request.build_absolute_uri("/").strip("/"),
        "DEPLOYMENT_VERSION": sha,
        "DEPLOYMENT_SOURCE_URL": f"https://github.com/tlambert03/FPbase/tree/{sha}" if sha else "",
    }
