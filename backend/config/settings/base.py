"""
Base settings for FPbase project.

For more information on this file, see
https://docs.djangoproject.com/en/dev/topics/settings/

For the full list of settings and their values, see
https://docs.djangoproject.com/en/dev/ref/settings/
"""
from pathlib import Path

import environ

ROOT_DIR = Path(__file__).resolve(strict=True).parent.parent.parent
APPS_DIR = ROOT_DIR / "fpbase"

# Load operating system environment variables and then prepare to use them
env = environ.Env()

# .env file, should load only in development environment
READ_DOT_ENV_FILE = env.bool("DJANGO_READ_DOT_ENV_FILE", default=False)

if READ_DOT_ENV_FILE:
    # Operating System Environment variables have precedence over variables defined in the .env file,
    # that is to say variables from the .env files will only be used if not defined
    # as environment variables.
    env_file = str(ROOT_DIR / ".env")
    print(f"Loading : {env_file}")
    env.read_env(env_file)
    print("The .env file has been loaded. See base.py for more information")

# APP CONFIGURATION
# ------------------------------------------------------------------------------
DJANGO_APPS = [
    # Default Django apps:
    "django.contrib.auth",
    "django.contrib.contenttypes",
    "django.contrib.sessions",
    "django.contrib.sites",
    "django.contrib.messages",
    "django.contrib.staticfiles",
    "django.contrib.sitemaps",
    # Useful template tags:
    "django.contrib.humanize",
    "dal",
    "dal_select2",
    # Admin
    "django.contrib.admin",
    # DRF
    "rest_framework",
    "drf_spectacular",
]

THIRD_PARTY_APPS = [
    "crispy_forms",  # Form layouts
    # "crispy_bootstrap4",
    "allauth",  # registration
    "allauth.account",  # registration
    "allauth.socialaccount",  # registration
    "allauth.socialaccount.providers.google",
    "allauth.socialaccount.providers.twitter",
    # 'allauth.socialaccount.providers.facebook',
    #  'allauth.socialaccount.providers.linkedin_oauth2',
    "captcha",
    "django_filters",
    "reversion",
    "reversion_compare",
    "avatar",
    "mptt",
]

# Apps specific for this project go here.
LOCAL_APPS = [
    # custom users app
    "fpbase.users.apps.UsersConfig",
    # Your stuff: custom apps go here
    "proteins.apps.ProteinsConfig",
    "references.apps.ReferencesConfig",
    "favit",
]

# See: https://docs.djangoproject.com/en/dev/ref/settings/#installed-apps
INSTALLED_APPS = DJANGO_APPS + THIRD_PARTY_APPS + LOCAL_APPS

# MIDDLEWARE CONFIGURATION
# ------------------------------------------------------------------------------
MIDDLEWARE = [
    "django.middleware.security.SecurityMiddleware",
    "django.contrib.sessions.middleware.SessionMiddleware",
    "django.middleware.common.CommonMiddleware",
    "django.middleware.csrf.CsrfViewMiddleware",
    "django.contrib.auth.middleware.AuthenticationMiddleware",
    "django.contrib.messages.middleware.MessageMiddleware",
    "django.middleware.clickjacking.XFrameOptionsMiddleware",
    "fpbase.middleware.CanonicalDomainMiddleware",
]

# MIGRATIONS CONFIGURATION
# ------------------------------------------------------------------------------
MIGRATION_MODULES = {"sites": "fpbase.contrib.sites.migrations"}

# DEBUG
# ------------------------------------------------------------------------------
# See: https://docs.djangoproject.com/en/dev/ref/settings/#debug
DEBUG = env.bool("DJANGO_DEBUG", False)

# FIXTURE CONFIGURATION
# ------------------------------------------------------------------------------
# See: https://docs.djangoproject.com/en/dev/ref/settings/#std:setting-FIXTURE_DIRS
FIXTURE_DIRS = (str(APPS_DIR / "fixtures"),)

# EMAIL CONFIGURATION
# ------------------------------------------------------------------------------
EMAIL_BACKEND = env("DJANGO_EMAIL_BACKEND", default="django.core.mail.backends.smtp.EmailBackend")

# MANAGER CONFIGURATION
# ------------------------------------------------------------------------------
# See: https://docs.djangoproject.com/en/dev/ref/settings/#admins
ADMINS = [("Talley Lambert", "talley.lambert+fpbase@gmail.com")]

# See: https://docs.djangoproject.com/en/dev/ref/settings/#managers
MANAGERS = [*ADMINS, ("Anna Jost", "anna_jost@hms.harvard.edu")]

EMAIL_SUBJECT_PREFIX = "[FPbase] "
SERVER_EMAIL = "FPbase <info@mg.fpbase.org>"
DEFAULT_FROM_EMAIL = "FPbase <info@mg.fpbase.org>"

# DATABASE CONFIGURATION
# ------------------------------------------------------------------------------
# See: https://docs.djangoproject.com/en/dev/ref/settings/#databases
# Uses django-environ to accept uri format
# See: https://django-environ.readthedocs.io/en/latest/#supported-types
DATABASES = {"default": env.db("DATABASE_URL", default="postgres:///fpbase")}
DATABASES["default"]["ATOMIC_REQUESTS"] = True

# https://docs.djangoproject.com/en/stable/ref/settings/#std:setting-DEFAULT_AUTO_FIELD
DEFAULT_AUTO_FIELD = "django.db.models.BigAutoField"

# GENERAL CONFIGURATION
# ------------------------------------------------------------------------------
# Local time zone for this installation. Choices can be found here:
# http://en.wikipedia.org/wiki/List_of_tz_zones_by_name
# although not all choices may be available on all operating systems.
# In a Windows environment this must be set to your system time zone.
TIME_ZONE = "UTC"

# See: https://docs.djangoproject.com/en/dev/ref/settings/#language-code
LANGUAGE_CODE = "en-us"

# See: https://docs.djangoproject.com/en/dev/ref/settings/#use-i18n
USE_I18N = True


# See: https://docs.djangoproject.com/en/dev/ref/settings/#use-tz
USE_TZ = True

# TEMPLATE CONFIGURATION
# ------------------------------------------------------------------------------
# See: https://docs.djangoproject.com/en/dev/ref/settings/#templates
TEMPLATES = [
    {
        # See: https://docs.djangoproject.com/en/dev/ref/settings/#std:setting-TEMPLATES-BACKEND
        "BACKEND": "django.template.backends.django.DjangoTemplates",
        # See: https://docs.djangoproject.com/en/dev/ref/settings/#template-dirs
        "DIRS": [str(APPS_DIR / "templates")],
        "OPTIONS": {
            # See: https://docs.djangoproject.com/en/dev/ref/settings/#template-debug
            "debug": DEBUG,
            # See: https://docs.djangoproject.com/en/dev/ref/settings/#template-loaders
            # https://docs.djangoproject.com/en/dev/ref/templates/api/#loader-types
            "loaders": [
                "django.template.loaders.filesystem.Loader",
                "django.template.loaders.app_directories.Loader",
            ],
            # See: https://docs.djangoproject.com/en/dev/ref/settings/#template-context-processors
            "context_processors": [
                "django.template.context_processors.debug",
                "django.template.context_processors.request",
                "django.contrib.auth.context_processors.auth",
                "django.template.context_processors.i18n",
                "django.template.context_processors.media",
                "django.template.context_processors.static",
                "django.template.context_processors.tz",
                "django.contrib.messages.context_processors.messages",
                # Your stuff: custom template context processors go here
                "fpbase.context_processors.api_keys",
                "fpbase.context_processors.canonical",
            ],
        },
    }
]

# See: http://django-crispy-forms.readthedocs.io/en/latest/install.html#template-packs
CRISPY_ALLOWED_TEMPLATE_PACKS = "bootstrap4"
CRISPY_TEMPLATE_PACK = "bootstrap4"


# STATIC FILE CONFIGURATION
# ------------------------------------------------------------------------------
# See: https://docs.djangoproject.com/en/dev/ref/settings/#static-root
STATIC_ROOT = str(ROOT_DIR / "staticfiles")

# See: https://docs.djangoproject.com/en/dev/ref/settings/#static-url
STATIC_URL = "/static/"

# See: https://docs.djangoproject.com/en/dev/ref/contrib/staticfiles/#std:setting-STATICFILES_DIRS
STATICFILES_DIRS = [str(ROOT_DIR.parent / "frontend" / "dist")]

# See: https://docs.djangoproject.com/en/dev/ref/contrib/staticfiles/#staticfiles-finders
STATICFILES_FINDERS = [
    "django.contrib.staticfiles.finders.FileSystemFinder",
    "django.contrib.staticfiles.finders.AppDirectoriesFinder",
]

INSTALLED_APPS.append("webpack_loader")

WEBPACK_LOADER = {
    "DEFAULT": {
        "CACHE": not DEBUG,
        "BUNDLE_DIR_NAME": "/",
        "STATS_FILE": str(ROOT_DIR.parent / "frontend" / "dist" / "webpack-stats.json"),
        "POLL_INTERVAL": 0.1,
        "TIMEOUT": None,
        "IGNORE": [r".*\.hot-update.js", r".+\.map"],
    }
}


# MEDIA CONFIGURATION
# ------------------------------------------------------------------------------
# See: https://docs.djangoproject.com/en/dev/ref/settings/#media-root
MEDIA_ROOT = str(APPS_DIR / "media")

# See: https://docs.djangoproject.com/en/dev/ref/settings/#media-url
MEDIA_URL = "/media/"

# URL Configuration
# ------------------------------------------------------------------------------
ROOT_URLCONF = "config.urls"

# See: https://docs.djangoproject.com/en/dev/ref/settings/#wsgi-application
WSGI_APPLICATION = "config.wsgi.application"

# PASSWORD STORAGE SETTINGS
# ------------------------------------------------------------------------------
# See https://docs.djangoproject.com/en/dev/topics/auth/passwords/#using-argon2-with-django
PASSWORD_HASHERS = [
    "django.contrib.auth.hashers.Argon2PasswordHasher",
    "django.contrib.auth.hashers.PBKDF2PasswordHasher",
    "django.contrib.auth.hashers.PBKDF2SHA1PasswordHasher",
    "django.contrib.auth.hashers.BCryptSHA256PasswordHasher",
    "django.contrib.auth.hashers.BCryptPasswordHasher",
]

# PASSWORD VALIDATION
# https://docs.djangoproject.com/en/dev/ref/settings/#auth-password-validators
# ------------------------------------------------------------------------------

AUTH_PASSWORD_VALIDATORS = [
    {"NAME": "django.contrib.auth.password_validation.UserAttributeSimilarityValidator"},
    {"NAME": "django.contrib.auth.password_validation.MinimumLengthValidator"},
    {"NAME": "django.contrib.auth.password_validation.CommonPasswordValidator"},
    {"NAME": "django.contrib.auth.password_validation.NumericPasswordValidator"},
]

# AUTHENTICATION CONFIGURATION
# ------------------------------------------------------------------------------
AUTHENTICATION_BACKENDS = [
    "django.contrib.auth.backends.ModelBackend",
    "allauth.account.auth_backends.AuthenticationBackend",
]

# Some really nice defaults
ACCOUNT_AUTHENTICATION_METHOD = "username"
ACCOUNT_EMAIL_REQUIRED = True
ACCOUNT_EMAIL_VERIFICATION = "mandatory"

ACCOUNT_ALLOW_REGISTRATION = env.bool("DJANGO_ACCOUNT_ALLOW_REGISTRATION", True)
ACCOUNT_ADAPTER = "fpbase.users.adapters.AccountAdapter"
SOCIALACCOUNT_ADAPTER = "fpbase.users.adapters.SocialAccountAdapter"
SOCIALACCOUNT_AUTO_SIGNUP = False

ACCOUNT_FORMS = {"signup": "fpbase.forms.CustomSignupForm"}

# Custom user app defaults
# Select the correct user model
AUTH_USER_MODEL = "users.User"
LOGIN_REDIRECT_URL = "users:redirect"
LOGIN_URL = "account_login"

# SLUGLIFIER
AUTOSLUG_SLUGIFY_FUNCTION = "slugify.slugify"

# django-rest-framework
# -------------------------------------------------------------------------------
# django-rest-framework - https://www.django-rest-framework.org/api-guide/settings/
REST_FRAMEWORK = {
    "DEFAULT_AUTHENTICATION_CLASSES": (
        "rest_framework.authentication.SessionAuthentication",
        "rest_framework.authentication.TokenAuthentication",
    ),
    "DEFAULT_PERMISSION_CLASSES": ("rest_framework.permissions.IsAuthenticated",),
    "DEFAULT_SCHEMA_CLASS": "drf_spectacular.openapi.AutoSchema",
}

# By Default swagger ui is available only to admin user(s). You can change permission classes to change that
# See more configuration options at https://drf-spectacular.readthedocs.io/en/latest/settings.html#settings
SPECTACULAR_SETTINGS = {
    "TITLE": "fpbase API",
    "DESCRIPTION": "Documentation of API endpoints of fpbase",
    "VERSION": "1.0.0",
    "SERVE_PERMISSIONS": ["rest_framework.permissions.IsAdminUser"],
}

# django-compressor
# ------------------------------------------------------------------------------
# INSTALLED_APPS += ['compressor']
# STATICFILES_FINDERS += ['compressor.finders.CompressorFinder']

# Location of root django.contrib.admin URL, use {% url 'admin:index' %}
ADMIN_URL = r"^admin/"

# See: https://docs.djangoproject.com/en/dev/ref/settings/#site-id
SITE_ID = 1
# CANONICAL_URL = env('CANONICAL_URL', default='https://www.fpbase.org')
CANONICAL_URL = env("CANONICAL_URL", default=None)


# AVATAR CONFIGURATION
# ------------------------------------------------------------------------------

AVATAR_AUTO_GENERATE_SIZES = (80, 36)
AVATAR_CACHE_ENABLED = True
AVATAR_GRAVATAR_DEFAULT = "identicon"
AVATAR_CLEANUP_DELETED = True
AVATAR_MAX_AVATARS_PER_USER = 8

# Your common stuff: Below this line define 3rd party library settings
# ------------------------------------------------------------------------------

MODERATION_MODERATORS = ("talley.lambert+fpbase@gmail.com",)

# v3 API for django-recaptcha
RECAPTCHA_PUBLIC_KEY = env("RECAPTCHA_V3_PUBLIC_KEY", default="")
RECAPTCHA_PRIVATE_KEY = env("RECAPTCHA_V3_PRIVATE_KEY", default="")
# NOCAPTCHA = True

GOOGLE_API_PRIVATE_KEY = env("GOOGLE_API_PRIVATE_KEY", default="").replace("#", "\n")
GOOGLE_API_CLIENT_EMAIL = env("GOOGLE_API_CLIENT_EMAIL", default=None)
GOOGLE_API_PRIVATE_KEY_ID = env("GOOGLE_API_PRIVATE_KEY_ID", default=None)

ALGOLIA_SUFFIX = "dev" if (DEBUG or ("staging" in env("SENTRY_PROJECT", default=""))) else "prod"
ALGOLIA_PUBLIC_KEY = "421b453d4f93e332ebd0c7f3ace29476"
ALGOLIA = {
    "APPLICATION_ID": "9WAWQMVNTB",
    "API_KEY": env("ALGOLIA_API_KEY", default=""),
    "INDEX_SUFFIX": ALGOLIA_SUFFIX,
}

if ALGOLIA["API_KEY"]:
    INSTALLED_APPS += ["algoliasearch_django"]

# CELERY_BROKER_URL = env('REDIS_URL', default='redis://localhost/')
CELERY_BROKER_URL = env("CLOUDAMQP_URL", default="amqp://localhost")
CELERY_RESULT_BACKEND = env("REDIS_URL", default="redis://localhost/")


INSTALLED_APPS += ["graphene_django"]
GRAPHENE = {"SCHEMA": "fpbase.schema.schema"}

# CORS
# -------

MIDDLEWARE = ["corsheaders.middleware.CorsMiddleware", *MIDDLEWARE]
INSTALLED_APPS += ["corsheaders"]
CORS_ORIGIN_WHITELIST = [
    "http://localhost:8000",
    "http://localhost:8080",
    "http://localhost:3000",
]
