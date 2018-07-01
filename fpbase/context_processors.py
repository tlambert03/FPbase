from django.conf import settings


def tracking_ids(request):
    return {
        'ga_tracking_id': getattr(settings, 'GA_TRACKING_ID', ''),
        'sentry_js_dsn': getattr(settings, 'SENTRY_JS_DSN', '')
    }


def canonical(request):
    return {'CANONICAL_URL': settings.CANONICAL_URL}
