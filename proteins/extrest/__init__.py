from ..models import Organism, Protein
from .entrez import get_gb_info
from .uniprot import get_uniprot_info


class ChangeSet:
    """Set of changes to make to a protein instance"""

    def __init__(self, obj):
        self.obj = obj
        self.id = obj.id
        self.changes = {}

    def __str__(self):
        s = f"{self.obj!s} ChangeSet:\n"
        for k, v in self.changes.items():
            if getattr(self.obj, k):
                s += f"\tCHANGE: {k} from {getattr(self.obj, k)} -> {v}\n"
            else:
                s += f"\tSET: {k} -> {v}\n"
        return s

    def __repr__(self):
        return f"<{self.obj} ChangeSet>"

    def __bool__(self):
        return bool(self.changes)

    def __add__(self, change):
        # where change is a tuple (attr, newval)
        if not (isinstance(change, tuple | list) and len(change) == 2):
            raise NotImplementedError("ChangeSet add epects a 2-tuple")

        attr, value = change
        if not hasattr(self.obj, attr):
            raise ValueError(f"Object {self.obj} does not have attribute {attr}")
        if attr in self.changes:
            if not value == self.changes[attr]:
                raise ValueError(f"Changeset received conflicting changes for field {attr}")
        else:
            self.changes[attr] = value
        return self

    def __iadd__(self, change):
        return self.__add__(change)

    def execute(self):
        if not self.changes:
            return
        p = Protein.objects.get(id=self.id)  # in case it's changed
        for attr, newval in self.changes.items():
            if attr == "parent_organism" and isinstance(newval, int):
                org, created = Organism.objects.get_or_create(id=newval)
                p.parent_organism = org
            else:
                setattr(p, attr, newval)
        p.save()


def compare_info(gb_info, up_info):
    mismatch = []
    if gb_info.get("gb_prot", False) and up_info.get("genbank", False):
        if gb_info["gb_prot"] not in up_info["genbank"]:
            mismatch.append("genbank mismatch gb_info: {}, uniprot: {}".format(gb_info["gb_prot"], up_info["genbank"]))
    if gb_info.get("uniprots", False) and up_info.get("uniprot", False):
        if up_info["uniprot"] not in gb_info["uniprots"]:
            mismatch.append(
                "uniprot mismatch uniprot: {}, gb_info: {}".format(up_info["uniprot"], gb_info["uniprots"])
            )
    if gb_info.get("seq", False) and up_info.get("seq", False):
        if not gb_info["seq"] == up_info["seq"]:
            mismatch.append("sequence mismatch gb_info:\n{}\n\nuniprot:\n{}".format(gb_info["seq"], up_info["seq"]))
    return mismatch


def integrity_check(protein):
    changes = ChangeSet(protein)
    warnings = []
    gbD = get_gb_info(protein.genbank) if protein.genbank else None
    upD = get_uniprot_info(protein.uniprot) if protein.uniprot else None
    if gbD and upD:
        [warnings.append(item) for item in compare_info(gbD, upD)]
    if gbD:
        # first make sure we're using the Protein accession (not the nuc)
        if gbD.get("gb_prot", False) and protein.genbank != gbD["gb_prot"]:
            changes += ("genbank", gbD["gb_prot"])
        # only add organism if non-existent
        if gbD.get("organism", False) and not protein.parent_organism:
            if not gbD["organism"] == 32630:
                changes += ("parent_organism", gbD["organism"])
        if gbD.get("ipg_id", False) and protein.ipg_id != gbD["ipg_id"]:
            changes += ("ipg_id", gbD["ipg_id"])
        if gbD.get("seq", False):
            if not protein.seq:
                changes += ("seq", gbD["seq"])
            else:
                if protein.seq not in gbD["seq"].upper():
                    warnings.append(
                        "{} sequence not found in genbank sequence!\n\tcurrent:{}\n\n\tgenbank:{}".format(
                            protein.name, protein.seq, gbD["seq"].upper()
                        )
                    )
        if len(gbD.get("uniprots", [])):
            if protein.uniprot:
                if protein.uniprot not in gbD["uniprots"]:
                    warnings.append(
                        "current uniprot id {} does not match the one from genbank {}".format(
                            protein.uniprot, gbD["uniprots"][0]
                        )
                    )
            else:
                changes += ("uniprot", gbD["uniprots"][0])
    if upD:
        if len(upD.get("genbank", [])) and not protein.genbank:
            changes += ("genbank", upD["genbank"][0])
        if upD.get("seq", False):
            if not protein.seq:
                changes += ("seq", upD["seq"])
            else:
                if protein.seq not in upD["seq"].upper():
                    warnings.append(
                        "{} sequence not found in genbank sequence!\n\tcurrent:{}\n\n\tgenbank:{}".format(
                            protein.name, protein.seq, upD["seq"].upper()
                        )
                    )

    return changes, warnings
