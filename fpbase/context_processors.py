from django.conf import settings


def api_keys(request):
    return {
        'ga_tracking_id': getattr(settings, 'GA_TRACKING_ID', ''),
        'sentry_dsn': getattr(settings, 'SENTRY_DSN', ''),
        'algolia_suffix': getattr(settings, 'ALGOLIA_SUFFIX', 'prod'),
        'algolia_public_key': getattr(settings, 'ALGOLIA_PUBLIC_KEY', ''),
        'algolia_app_id': getattr(settings, 'ALGOLIA', {}).get('APPLICATION_ID')
    }


def canonical(request):
    return {'CANONICAL_URL': settings.CANONICAL_URL}
