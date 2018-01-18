from django.conf import settings


def ga_tracking_id(request):
    return {'ga_tracking_id': getattr(settings, 'GA_TRACKING_ID', '')}
