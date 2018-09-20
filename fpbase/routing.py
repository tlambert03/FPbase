from channels.routing import ProtocolTypeRouter, URLRouter
from channels.auth import AuthMiddlewareStack
import proteins.routing

application = ProtocolTypeRouter({
    # Empty for now (http->django views is added by default)
    # (http->django views is added by default)
    'websocket': AuthMiddlewareStack(
        URLRouter(
            proteins.routing.websocket_urlpatterns
        )
    ),
})
