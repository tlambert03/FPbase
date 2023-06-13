from allauth.account.signals import user_signed_up
from django.core.mail import mail_admins
from django.dispatch import receiver


@receiver(user_signed_up)
def notify_admins_user_signed_up(sender, **kwargs):
    if kwargs.get("user", False):
        mail_admins(
            "New user created: {}".format(kwargs.get("user")),
            "User: {} created a new account".format(kwargs.get("user")),
            fail_silently=True,
        )
