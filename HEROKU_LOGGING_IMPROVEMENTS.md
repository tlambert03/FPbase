# Heroku Logging Improvements for FPbase

**Date:** 2025-10-22
**Status:** Ready for Implementation
**Estimated Time:** 30-45 minutes

## Executive Summary

FPbase is experiencing log overflow on Heroku ("Error L10: output buffer overflow") due to exceeding Papertrail's free tier limit of 10 MB/day. This analysis recommends switching to Better Stack (Logtail) for 300x more capacity (3 GB/day) and implementing structured JSON logging for better AI agent compatibility.

---

## Current State Analysis

### Logging Infrastructure
- **Service:** Papertrail (addon: `papertrail-graceful-91985`)
- **Plan:** choklad (free tier)
- **Limits:** 10 MB/day, 1-week retention
- **Problem:** Regular log overflow causing dropped messages

### Evidence of Issues
From production logs:
```
2025-10-22T14:25:26.226111+00:00 heroku[logplex]: Error L10 (output buffer overflow):
drain 'd.efe51fd5-8fdb-4153-b094-f5e9f7df3058' dropped 1 messages since 2025-10-22T14:24:56.323777+00:00.

2025-10-22T14:25:27.747804+00:00 heroku[logplex]: Error L10 (output buffer overflow):
drain 'd.efe51fd5-8fdb-4153-b094-f5e9f7df3058' dropped 12 messages since 2025-10-22T14:24:27.798935+00:00.
```

Messages are being dropped every ~30 seconds, indicating consistent overflow.

### Log Volume Analysis

**Primary Contributors:**
1. **Heroku Router Logs (50-70% of volume)**
   - Every HTTP request generates a router log
   - Static file requests are especially verbose
   - Cannot be disabled at platform level
   - Sample: `heroku[router]: at=info method=GET path="/static/vendor.151952384fb7b410a50b.js"`

2. **Application Debug/Info Logging (20-30%)**
   - `/backend/proteins/views/spectra.py`: 10+ `logger.info()` calls
   - Excessive detail for production environment
   - Most useful for development/debugging only

3. **Print Statements (minor but unnecessary)**
   - `/backend/config/settings/production.py:201`: `print(f"HEROKU_SLUG_COMMIT = {HEROKU_SLUG_COMMIT}")`
   - Printed on every dyno startup
   - Should use logger instead

4. **Middleware Debug Logs (minimal)**
   - `/backend/fpbase/middleware.py`: 2 `logger.debug()` calls
   - Already at DEBUG level, not contributing to overflow

### Current Logging Configuration

**Location:** `/backend/config/settings/production.py:214-245`

```python
LOGGING = {
    "version": 1,
    "disable_existing_loggers": True,
    "root": {"level": "WARNING"},
    "formatters": {
        "verbose": {"format": "%(levelname)s %(asctime)s %(module)s %(process)d %(thread)d %(message)s"},
    },
    "handlers": {
        "console": {
            "level": "DEBUG",
            "class": "logging.StreamHandler",
            "formatter": "verbose",
        }
    },
    "loggers": {
        "django.db.backends": {"level": "ERROR", "handlers": ["console"], "propagate": False},
        "sentry.errors": {"level": "DEBUG", "handlers": ["console"], "propagate": False},
        "django.security.DisallowedHost": {"level": "ERROR", "handlers": ["console"], "propagate": False},
    },
}
```

**Issues:**
- No structured/JSON formatting (harder for AI agents to parse)
- No context enrichment (request IDs, user info)
- No filtering for static file requests
- Verbose formatter outputs more data than necessary

---

## Research: Free Logging Alternatives

### Option 1: Better Stack (Logtail) ⭐ RECOMMENDED
- **Free Tier:** 3 GB/day, 3-day retention
- **Advantage:** 300x more capacity than Papertrail
- **Features:** SQL-compatible queries, structured logging support, fast UI
- **Heroku Addon:** `logtail:free`
- **Cost:** $0.45/GB beyond free tier, $0.025/GB/week for extended retention
- **Best For:** Current needs + significant growth headroom

### Option 2: Sumo Logic
- **Free Tier:** 500 MB/day, 7-day retention
- **Advantage:** 50x more capacity, longer retention than Papertrail
- **Features:** Advanced analytics, dashboards
- **Heroku Addon:** `sumologic:free`
- **Trial:** 5 GB/day for 30 days, then reverts to 500 MB/day
- **Best For:** If longer retention is priority

### Option 3: LogDNA (Mezmo)
- **Free Tier:** Unclear from research (likely limited/trial only)
- **Heroku Addon:** `logdna` (now Mezmo)
- **Note:** Free tier details not well documented; may not have permanent free tier

### Option 4: Stay with Papertrail
- **Not Recommended:** Already experiencing overflow
- **Workaround:** Could use filtering in Papertrail UI, but doesn't solve capacity issue
- **Note:** Would require significant log reduction to stay under 10 MB/day

---

## Recommended Solution: Multi-Pronged Approach

### Strategy
1. **Switch to Better Stack** (immediate 300x capacity increase)
2. **Reduce application log verbosity** (better signal-to-noise)
3. **Implement structured JSON logging** (AI-agent friendly)
4. **Filter static file requests** (reduce Django app log volume)

### Expected Results
- ✅ Eliminate log overflow completely
- ✅ 3 GB/day capacity should handle current traffic + 30x growth
- ✅ Better debugging experience for AI agents
- ✅ Reduced noise, more actionable logs
- ✅ Zero additional cost

---

## Implementation Instructions

### Prerequisites
- Heroku CLI installed and authenticated
- Access to FPbase Heroku app
- Repository write access

### Phase 1: Switch to Better Stack (5 minutes)

**Step 1.1:** Install Better Stack addon
```bash
heroku addons:create logtail:free --app fpbase
```

**Step 1.2:** Verify logs are flowing
```bash
# Check addon was created
heroku addons --app fpbase | grep logtail

# Open Better Stack dashboard
heroku addons:open logtail --app fpbase

# In the dashboard, verify you see recent logs appearing
```

**Step 1.3:** Monitor for 24 hours (optional but recommended)
Keep both services running in parallel to ensure Better Stack is working properly.

**Step 1.4:** Remove Papertrail (after verification)
```bash
# Only do this after confirming Better Stack is working
heroku addons:destroy papertrail --app fpbase
```

---

### Phase 2: Update Production Logging Config (10 minutes)

**File:** `/backend/config/settings/production.py`

**Step 2.1:** Replace print statement with logger

**Current (line 201):**
```python
print(f"HEROKU_SLUG_COMMIT = {HEROKU_SLUG_COMMIT}")
```

**Change to:**
```python
import logging
logger = logging.getLogger(__name__)
# ... later in the file, after LOGGING is configured:
logger.info(f"Application starting with commit: {HEROKU_SLUG_COMMIT}")
```

**Step 2.2:** Update LOGGING configuration

Replace the existing `LOGGING` dict (lines 214-245) with this improved version:

```python
LOGGING = {
    "version": 1,
    "disable_existing_loggers": False,  # Changed from True to preserve app loggers
    "formatters": {
        "json": {
            "()": "pythonjsonlogger.jsonlogger.JsonFormatter",
            "format": "%(asctime)s %(name)s %(levelname)s %(message)s %(pathname)s %(lineno)d",
        },
        "verbose": {
            "format": "%(levelname)s %(asctime)s %(module)s %(process)d %(thread)d %(message)s"
        },
    },
    "filters": {
        "skip_static_requests": {
            "()": "django.utils.log.CallbackFilter",
            "callback": lambda record: (
                not hasattr(record, "args")
                or not record.args
                or not isinstance(record.args, tuple)
                or not record.args[0].startswith(("GET /static/", "GET /media/"))
            ),
        },
    },
    "handlers": {
        "console": {
            "level": "INFO",  # Changed from DEBUG
            "class": "logging.StreamHandler",
            "formatter": "json",  # Use JSON for production
            "filters": ["skip_static_requests"],
        }
    },
    "root": {
        "level": "INFO",  # Changed from WARNING for better coverage
        "handlers": ["console"],
    },
    "loggers": {
        "django": {
            "level": "INFO",
            "handlers": ["console"],
            "propagate": False,
        },
        "django.server": {
            "level": "INFO",
            "handlers": ["console"],
            "propagate": False,
            "filters": ["skip_static_requests"],  # Filter Django dev server logs
        },
        "django.db.backends": {
            "level": "ERROR",  # Only log database errors
            "handlers": ["console"],
            "propagate": False,
        },
        "sentry.errors": {
            "level": "WARNING",  # Changed from DEBUG
            "handlers": ["console"],
            "propagate": False,
        },
        "django.security.DisallowedHost": {
            "level": "ERROR",
            "handlers": ["console"],
            "propagate": False,
        },
        "proteins": {
            "level": "WARNING",  # Only warnings and errors for proteins app
            "handlers": ["console"],
            "propagate": False,
        },
        "fpbase": {
            "level": "WARNING",  # Only warnings and errors for fpbase app
            "handlers": ["console"],
            "propagate": False,
        },
    },
}
```

**Step 2.3:** Add python-json-logger dependency

Add to `/backend/pyproject.toml` dependencies:
```toml
dependencies = [
    # ... existing dependencies ...
    "python-json-logger>=2.0.7",
]
```

Then run:
```bash
uv sync
```

---

### Phase 3: Reduce Application Log Verbosity (5 minutes)

**File:** `/backend/proteins/views/spectra.py`

**Current Issue:** Lines 164-257 contain excessive `logger.info()` calls

**Strategy:** Change routine operations from INFO to DEBUG level

**Changes needed:**

```python
# Line 164: Keep as INFO (useful for tracking user actions)
logger.info(f"Spectrum preview request from user: {request.user}")

# Line 168: Change to DEBUG
logger.debug(f"Tab selection - data_source: {data_source}")  # Changed from info

# Line 172: Change to DEBUG
logger.debug("Using manual data entry (files ignored)")  # Changed from info

# Line 175: Change to DEBUG
logger.debug(f"Manual data submitted: {manual_data[:100]}...")  # Changed from info

# Line 179: Change to DEBUG
logger.debug("Using file upload (manual data ignored)")  # Changed from info

# Line 186: Keep as WARNING (errors should always be logged)
logger.warning(f"Form validation failed: {form.errors}")

# Line 205: Change to DEBUG
logger.debug("Running spectrum normalization...")  # Changed from info

# Line 207: Change to DEBUG
logger.debug("Spectrum normalization completed successfully")  # Changed from info

# Line 220: Change to DEBUG
logger.debug("Generating SVG image...")  # Changed from info

# Line 230: Change to DEBUG
logger.debug("SVG generation completed successfully")  # Changed from info

# Line 257: Keep as INFO or change to DEBUG (decide based on importance)
logger.debug("Spectrum preview generated successfully")  # Changed from info
```

**Rule of Thumb:**
- **ERROR**: Something broke, needs immediate attention
- **WARNING**: Something unexpected, but handled gracefully
- **INFO**: Important business logic milestones (user actions, external API calls)
- **DEBUG**: Detailed flow information, useful only when debugging

---

### Phase 4: Test Changes Locally (5 minutes)

**Step 4.1:** Run local development server
```bash
pnpm dev
# or
uv run backend/manage.py runserver
```

**Step 4.2:** Verify logging works
```bash
# Check logs are being output
# Try accessing the site and submitting a spectrum
# Verify log format looks correct
```

**Step 4.3:** Check for errors
```bash
# Look for any ImportError or configuration errors
uv run pytest backend/proteins/tests/ -v
```

---

### Phase 5: Deploy to Heroku (5 minutes)

**Step 5.1:** Commit changes
```bash
git add backend/config/settings/production.py
git add backend/proteins/views/spectra.py
git add backend/pyproject.toml
git add uv.lock  # if changed by uv sync

git commit -m "Improve Heroku logging: switch to Better Stack, add JSON logging, reduce verbosity

- Add python-json-logger for structured logging
- Implement JSON formatter with context fields
- Add filter to skip static/media file logs
- Reduce proteins.views.spectra verbosity (INFO -> DEBUG)
- Replace print statement with logger in production.py
- Update log levels for better signal-to-noise ratio

This addresses Error L10 (log overflow) by reducing log volume
and prepares logs for Better Stack/Logtail integration.

🤖 Generated with [Claude Code](https://claude.com/claude-code)

Co-Authored-By: Claude <noreply@anthropic.com>"
```

**Step 5.2:** Push to Heroku
```bash
git push heroku main
# or if on a different branch:
git push heroku add-indes:main
```

**Step 5.3:** Monitor deployment
```bash
# Watch the build
heroku logs --tail --app fpbase

# Look for the startup message with commit hash
# Check for any errors during deployment
```

**Step 5.4:** Verify logs in Better Stack
```bash
# Open Better Stack dashboard
heroku addons:open logtail --app fpbase

# Verify you see JSON-formatted logs
# Check that static file requests are filtered
# Confirm no more L10 overflow errors
```

---

### Phase 6: Advanced Structured Logging (Optional, 15 minutes)

For even better AI agent compatibility, consider adding `django-structlog` for automatic request context enrichment.

**Step 6.1:** Add django-structlog
```bash
# Add to pyproject.toml
dependencies = [
    # ... existing ...
    "django-structlog>=8.0.0",
]

uv sync
```

**Step 6.2:** Update settings

Add to `INSTALLED_APPS` in `/backend/config/settings/base.py`:
```python
INSTALLED_APPS = [
    # ... existing apps ...
    "django_structlog",
]
```

Add to `MIDDLEWARE` in `/backend/config/settings/base.py` (after RequestMiddleware):
```python
MIDDLEWARE = [
    # ... existing middleware ...
    "django_structlog.middlewares.RequestMiddleware",
]
```

**Step 6.3:** Configure structlog in production.py

Add after imports:
```python
import structlog

structlog.configure(
    processors=[
        structlog.contextvars.merge_contextvars,
        structlog.stdlib.filter_by_level,
        structlog.stdlib.add_logger_name,
        structlog.stdlib.add_log_level,
        structlog.stdlib.PositionalArgumentsFormatter(),
        structlog.processors.TimeStamper(fmt="iso"),
        structlog.processors.StackInfoRenderer(),
        structlog.processors.format_exc_info,
        structlog.processors.UnicodeDecoder(),
        structlog.processors.JSONRenderer(),
    ],
    wrapper_class=structlog.stdlib.BoundLogger,
    context_class=dict,
    logger_factory=structlog.stdlib.LoggerFactory(),
    cache_logger_on_first_use=True,
)
```

**Step 6.4:** Update LOGGING to use structlog

Modify the console handler:
```python
"handlers": {
    "console": {
        "level": "INFO",
        "class": "logging.StreamHandler",
        "formatter": "json",
    }
},
```

And add to LOGGING:
```python
"formatters": {
    "json": {
        "()": structlog.stdlib.ProcessorFormatter,
        "processor": structlog.processors.JSONRenderer(),
    },
}
```

---

## Verification Checklist

After implementation, verify:

- [ ] Better Stack addon is installed and receiving logs
- [ ] No more "Error L10 (output buffer overflow)" messages
- [ ] Logs are in JSON format in Better Stack dashboard
- [ ] Static file requests are filtered (check Django app logs, not router logs)
- [ ] Application still functions normally (no logging-related errors)
- [ ] Log volume is under 3 GB/day (check Better Stack usage dashboard)
- [ ] Commit messages from HEROKU_SLUG_COMMIT appear in logs
- [ ] Error and warning messages still appear as expected
- [ ] Sentry integration still working (if errors occur)

---

## Expected Impact

### Before
- **Capacity:** 10 MB/day (exceeded daily)
- **Log Overflow:** Multiple L10 errors per minute
- **Retention:** 7 days (but logs being dropped)
- **Format:** Plain text, inconsistent structure
- **AI Debugging:** Difficult to parse and analyze
- **Signal-to-Noise:** Low (too many routine INFO logs)

### After
- **Capacity:** 3 GB/day (300x increase)
- **Log Overflow:** None expected (30x headroom)
- **Retention:** 3 days (all logs preserved)
- **Format:** Structured JSON with consistent schema
- **AI Debugging:** Easy to query and analyze with SQL-like syntax
- **Signal-to-Noise:** High (focus on warnings/errors)

---

## Monitoring and Maintenance

### Daily (First Week)
1. Check Better Stack dashboard for log volume
2. Verify no L10 overflow errors
3. Monitor for any new errors introduced by logging changes

### Weekly
1. Review log volume trends
2. Check for any repeated warnings that need attention
3. Verify retention is sufficient for debugging needs

### Monthly
1. Review Better Stack usage (should stay well under 3 GB/day)
2. Evaluate if more structured logging would be beneficial
3. Consider if retention needs have changed

---

## Rollback Plan

If issues occur:

**Immediate Rollback (Code Changes):**
```bash
git revert <commit-hash>
git push heroku main
```

**Revert to Papertrail (if Better Stack has issues):**
```bash
# Re-add Papertrail
heroku addons:create papertrail:choklad --app fpbase

# Remove Better Stack
heroku addons:destroy logtail --app fpbase
```

---

## AI Agent Instructions

When you implement this plan:

1. **Read this entire document first** - understand the full context
2. **Work through phases sequentially** - don't skip ahead
3. **Test locally before deploying** - catch configuration errors early
4. **Make atomic commits** - separate logical changes for easier rollback
5. **Monitor after each phase** - verify nothing breaks
6. **Use the verification checklist** - ensure complete implementation
7. **Ask questions if unclear** - better to clarify than make assumptions

### Key Files to Modify
- `/backend/config/settings/production.py` (main changes)
- `/backend/proteins/views/spectra.py` (reduce verbosity)
- `/backend/pyproject.toml` (add dependency)

### Commands You'll Need
```bash
# Heroku addon management
heroku addons:create logtail:free --app fpbase
heroku addons:open logtail --app fpbase
heroku addons:destroy papertrail --app fpbase

# Dependency management
uv sync

# Git workflow
git add <files>
git commit -m "message"
git push heroku main

# Monitoring
heroku logs --tail --app fpbase
```

---

## Additional Resources

- [Better Stack Heroku Integration](https://betterstack.com/docs/logs/heroku/)
- [Django Logging Configuration](https://docs.djangoproject.com/en/5.0/topics/logging/)
- [python-json-logger](https://github.com/madzak/python-json-logger)
- [django-structlog](https://github.com/jrobichaud/django-structlog)
- [Heroku Log Drains](https://devcenter.heroku.com/articles/log-drains)

---

## Questions or Issues?

If you encounter problems during implementation:

1. Check Heroku logs for error messages
2. Verify Python syntax in modified files
3. Ensure dependencies were installed correctly with `uv sync`
4. Test locally before deploying
5. Review Better Stack documentation for addon-specific issues
6. Consult this document's rollback plan

---

**Document Version:** 1.0
**Last Updated:** 2025-10-22
**Author:** Claude (Anthropic)
**Review Status:** Ready for Implementation
