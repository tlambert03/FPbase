from django.conf import settings


def ga_tracking_id(request):
    return {'ga_tracking_id': getattr(settings, 'GA_TRACKING_ID', '')}


def canonical(request):
    return {'CANONICAL_URL': settings.CANONICAL_URL}
