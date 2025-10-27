# End-to-End Tests

This directory contains end-to-end (e2e) tests for FPbase using Playwright and Django's `live_server` fixture.

## Running E2E Tests

E2E tests are **excluded from the default test run** to prevent database cleanup conflicts with unit tests.

```bash
# Run ONLY e2e tests
uv run pytest backend/tests_e2e/

# Run a specific e2e test
uv run pytest backend/tests_e2e/test_e2e.py::test_main_page

# Run with headed browser (to see what's happening)
uv run pytest backend/tests_e2e/ --headed

# Run with specific browser
uv run pytest backend/tests_e2e/ --browser chromium
uv run pytest backend/tests_e2e/ --browser firefox
```

## Why E2E Tests Are Separate

### Database Cleanup Differences

**Unit Tests:**
- Use Django's `TestCase` which wraps each test in a transaction
- After each test: **ROLLBACK** transaction (fast, reliable)
- Safe to run many tests in parallel

**E2E Tests (this directory):**
- Use `live_server` which runs Django in a **separate thread**
- Transactions can't span threads, so rollback is impossible
- After each test: **TRUNCATE all database tables** (slower)
- Must be isolated from unit tests to avoid conflicts

### Key Constraints

1. **No session-scoped database fixtures**: Data created in session-scoped fixtures will be wiped after the first test
2. **Each test must be independent**: Don't rely on data from previous tests
3. **Use factories for test data**: Create data in each test function or in function-scoped fixtures

## Configuration

All e2e-specific configuration is in `conftest.py`:

- **`pytestmark`**: Marks all tests with `django_db(transaction=True)` (required for live_server)
- **`_build_frontend_assets`**: Module-scoped fixture to build webpack assets once
- **`page`**: Configures Playwright page with sensible defaults (timeout, viewport)
- **`DJANGO_ALLOW_ASYNC_UNSAFE`**: Required for Playwright's async event loop to work with Django ORM

## References

- [pytest-playwright-and-django blog post](https://blog.tmk.name/2025/04/06/pytest-playwright-and-django/)
- [pytest-django database documentation](https://pytest-django.readthedocs.io/en/latest/database.html)
- [Playwright Python documentation](https://playwright.dev/python/)

## Writing New E2E Tests

```python
def test_my_feature(live_server: LiveServer, page: Page) -> None:
    """Test description."""
    # Navigate to page
    page.goto(f"{live_server.url}/my-feature/")

    # Interact with page
    page.locator("button").click()

    # Assert on results
    assert page.locator("h1").text_content() == "Expected Title"
```

### Best Practices

1. **Use proper waits**: Don't use `time.sleep()`, use Playwright's built-in waits
   ```python
   # Good
   page.wait_for_selector(".result")
   page.locator("button").wait_for(state="visible")

   # Bad
   import time
   time.sleep(2)
   ```

2. **Create test data in the test**: Use factories to create exactly what you need
   ```python
   from proteins.factories import ProteinFactory

   def test_protein_page(live_server, page):
       protein = ProteinFactory(name="TestProtein")
       page.goto(f"{live_server.url}{protein.get_absolute_url()}")
       # ...
   ```

3. **Clean up after yourself** (if needed): Test isolation is automatic via table truncation, but if you create external resources (files, etc), clean them up in a fixture

4. **Use descriptive locators**: Prefer semantic selectors over CSS selectors
   ```python
   # Good
   page.get_by_role("button", name="Submit")
   page.get_by_label("Username")
   page.get_by_text("Welcome")

   # Acceptable
   page.locator("#submit-btn")

   # Less ideal (brittle)
   page.locator("div > div > button.btn.btn-primary")
   ```
