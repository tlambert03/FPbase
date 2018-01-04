from moderation import moderation
from moderation.moderator import GenericModerator
from .models import Protein, State


class MyBaseModerator(GenericModerator):
    notify_user = True
    notify_moderator = True
    auto_approve_for_superusers = True
    fields_exclude = ['created', 'modified']


class ProteinModerator(MyBaseModerator):
    pass


class StateModerator(MyBaseModerator):
    pass


moderation.register(Protein, ProteinModerator)  # Uses default moderation settings
moderation.register(State, StateModerator)  # Uses default moderation settings
