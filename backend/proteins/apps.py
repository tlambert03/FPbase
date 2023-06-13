from django.apps import AppConfig


class ProteinsConfig(AppConfig):
    name = "proteins"

    def ready(self):
        # Makes sure all signal handlers are connected
        from . import handlers  # noqa
