import factory.declarations


class UserFactory(factory.django.DjangoModelFactory):
    username = factory.declarations.Sequence(lambda n: f"user-{n}")
    email = factory.declarations.Sequence(lambda n: f"user-{n}@example.com")
    password = factory.declarations.PostGenerationMethodCall("set_password", "password")

    class Meta:
        model = "users.User"
        django_get_or_create = ("username",)
