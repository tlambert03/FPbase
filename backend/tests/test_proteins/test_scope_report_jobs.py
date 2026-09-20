"""Tests for the job-control actions on the microscope report page."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest
from django.urls import reverse

from proteins.factories import MicroscopeFactory

AJAX = {"HTTP_X_REQUESTED_WITH": "XMLHttpRequest"}


@pytest.fixture
def celery_app():
    with patch("proteins.views.microscope.app") as app:
        app.control.inspect.return_value.active.return_value = None
        app.AsyncResult.return_value = MagicMock(ready=lambda: False, info=None)
        yield app


@pytest.fixture
def delay():
    with patch("proteins.views.microscope.calculate_scope_report.delay") as delay:
        delay.return_value.id = "real-job-id"
        yield delay


def _post(client, scope, **data):
    url = reverse("proteins:microscope-report", args=[scope.id])
    return client.post(url, data, **AJAX)


@pytest.mark.django_db
@pytest.mark.parametrize("action", ["cancel", "check"])
def test_arbitrary_job_id_is_rejected(client, celery_app, action):
    response = _post(client, MicroscopeFactory(), action=action, job_id="someone-elses-job")
    assert response.json() == {"status": 404}
    celery_app.AsyncResult.assert_not_called()


@pytest.mark.django_db
@pytest.mark.parametrize("action", ["cancel", "check"])
def test_job_token_only_valid_for_its_microscope(client, celery_app, delay, action):
    scope, other_scope = MicroscopeFactory(), MicroscopeFactory()
    token = _post(client, scope, action="update").json()["job"]
    delay.assert_called_once_with(scope.id, outdated_ids=None)
    assert token != "real-job-id"

    assert _post(client, other_scope, action=action, job_id=token).json() == {"status": 404}
    celery_app.AsyncResult.assert_not_called()

    response = _post(client, scope, action=action, job_id=token)
    assert response.json()["status"] == 200
    celery_app.AsyncResult.assert_called_once_with("real-job-id")
    if action == "cancel":
        celery_app.AsyncResult.return_value.revoke.assert_called_once_with(terminate=True)


@pytest.mark.django_db
def test_update_returns_token_for_already_running_job(client, celery_app, delay):
    scope = MicroscopeFactory()
    celery_app.control.inspect.return_value.active.return_value = {
        "worker": [
            {
                "id": "running-id",
                "name": "proteins.tasks.calculate_scope_report",
                "args": [scope.id],
            }
        ]
    }
    token = _post(client, scope, action="update").json()["job"]
    delay.assert_not_called()

    _post(client, scope, action="cancel", job_id=token)
    celery_app.AsyncResult.assert_called_once_with("running-id")


@pytest.mark.django_db
def test_report_post_unknown_microscope_404(client, celery_app, delay):
    url = reverse("proteins:microscope-report", args=["doesnotexist"])
    response = client.post(url, {"action": "update"}, **AJAX)
    assert response.status_code == 404
    delay.assert_not_called()
