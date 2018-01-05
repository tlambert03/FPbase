from moderation import moderation
from moderation.moderator import GenericModerator
from .models import Protein, State


class MyBaseModerator(GenericModerator):
    notify_user = False
    notify_moderator = False
    auto_approve_for_superusers = True
    auto_approve_for_staff = True
    auto_reject_for_anonymous = False
    fields_exclude = ['created', 'modified']


class ProteinModerator(MyBaseModerator):
    manager_names = ['objects']


class StateModerator(MyBaseModerator):
    pass


moderation.register(Protein, ProteinModerator)  # Uses default moderation settings
moderation.register(State, StateModerator)  # Uses default moderation settings
