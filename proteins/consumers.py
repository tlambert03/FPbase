from channels.generic.websocket import WebsocketConsumer
from django.core.cache import cache
import json
from .util.helpers import forster_list


class FRETConsumer(WebsocketConsumer):
    def connect(self):
        print("FRETConsumer accepted")
        self.accept()
        L = cache.get('forster_list')
        self.send(text_data=json.dumps({'message': L}))
        if not L:
            L = forster_list()
            cache.set('forster_list', L, 60 * 60 * 24 * 7)
            self.send(text_data=json.dumps({'message': L}))

    def disconnect(self, close_code):
        pass

    def receive(self, text_data):
        pass
