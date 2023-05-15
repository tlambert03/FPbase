import contextlib
from collections import OrderedDict, defaultdict

from reversion.models import Revision, Version
from reversion.revisions import transaction

from references.models import Reference


def lists_are_different(a, b):
    return (set(a) != set(b)) or (len(a) != len(b))


def listdiff(a, b):
    aset = set(a)
    bset = set(b)
    d = {}
    if aset.difference(bset):
        d["removed"] = ", ".join([str(x) for x in aset.difference(bset)])
    if bset.difference(aset):
        d["added"] = ", ".join([str(x) for x in bset.difference(aset)])
    if (aset == bset) and len(a) != len(b):
        d["changed"] = f"{a}->{b}"
    return d or None


def dictdiff(a, b, ignoreKeys=()):
    if not (isinstance(a, dict) and isinstance(b, dict)):
        if isinstance(a, list) and isinstance(b, list):
            return listdiff(a, b)
        # elif isinstance(a, (int, float, str, bool)) and isinstance(b, (int, float, str, bool)):
        #    return {'changed': (a, b)}
        elif a.__class__.__name__ == "FPSeq" or b.__class__.__name__ == "FPSeq":
            if a and b:
                return {"changed": str(a.mutations_to(b))}
            if a and not b:
                return {"removed": str(a)}
            if b and not a:
                return {"added": str(b)}
        elif isinstance(a, bool) and isinstance(b, bool):
            return {"changed": f"{a}->{b}"}
        else:
            if a and not b:
                return {"removed": a}
            if b and not a:
                return {"added": b}
            else:
                return {"changed": f"{a}->{b}"}

    res = {}
    for key in sorted(set(list(a.keys()) + list(b.keys()))):
        if key in ignoreKeys:
            continue
        v1 = a.get(key, None)
        v2 = b.get(key, None)
        if v1 and not v2:
            res[key] = {"removed": v1}
        if v2 and not v1:
            res[key] = {"added": v2}
        if v1 != v2:
            if isinstance(v1, list) and isinstance(v2, list):
                if lists_are_different(v1, v2):
                    res[key] = dictdiff(v1, v2, ignoreKeys)
            else:
                res[key] = dictdiff(v1, v2, ignoreKeys)
    return res


class UndoRollback(Exception):
    def __init__(self, object):
        self.object = object


def old_object(ver):
    try:
        with transaction.atomic(using=ver.db, savepoint=False):
            # Revert to this ver
            ver.revision.revert(delete=True)
            Model = ver.content_type.model_class()
            # get the old object
            old = Model.objects.prefetch_related("states").get(id=ver.object_id)
            # Raise the error to roll back the changes
            raise UndoRollback(old)
    except UndoRollback as ex:
        return ex.object


def get_history(obj, ignoreKeys=()):
    rev_ids = sorted(Version.objects.get_for_object(obj).values_list("revision__id", flat=True))

    revisions = list(
        Revision.objects.filter(id__in=rev_ids)
        .order_by("date_created")
        .prefetch_related("version_set", "user", "version_set__content_type")
    )

    object_versions = defaultdict(list)
    object_revisions = defaultdict(list)
    object_reprs = {}
    for rev in revisions:
        for v in rev.version_set.all():
            object_repr = (v.object_id, v.content_type_id)
            object_versions[object_repr].append(v)
            object_revisions[object_repr].append(rev)
            object_reprs[object_repr] = f"{v.content_type.name} {v.object_repr}"

    changes = OrderedDict([(rev, defaultdict(list)) for rev in revisions[1:]])
    # changes[revisions[0]] = {"initial revision"
    for object_repr, versions in object_versions.items():
        rep = object_reprs[object_repr]
        for n in range(len(versions) - 1):
            diffs = dictdiff(versions[n].field_dict, versions[n + 1].field_dict, ignoreKeys)
            if diffs:
                revision = versions[n + 1].revision
                for field, actions in diffs.items():
                    for action, result in actions.items():
                        if field == "references":
                            _result = []
                            for ref_id in result.split(", "):
                                with contextlib.suppress(Reference.DoesNotExist):
                                    ref = Reference.objects.get(id=ref_id)
                                    _result.append(ref.citation)
                            result = ", ".join(_result)
                        changes[revision][rep].append((action, field, result))
        revs = object_revisions[object_repr]
        if max(r.id for r in revs) < max(rev_ids):
            del_at = next(i for i in revisions if i.id > max(r.id for r in revs))
            changes[del_at][rep].append(("removed", None, None))
        if min(r.id for r in revs) > min(rev_ids):
            changes[revs[0]][rep].append(("added", None, None))

    for _c, v in changes.items():
        v.default_factory = None

    return changes
