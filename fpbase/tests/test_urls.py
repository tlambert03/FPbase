from django.urls import reverse
from django import test
from config.urls import urlpatterns


class UrlsTest(test.TestCase):

    def test_responses(self):
        for url in urlpatterns:
            if hasattr(url, 'name'):
                print(url)
                response = self.client.get(reverse(url.name))
                self.assertEqual(response.status_code, 200)
