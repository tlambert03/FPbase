from __future__ import annotations

from unittest.mock import patch

import pytest

from proteins.extrest import entrez


@pytest.mark.parametrize(("id_list", "expected"), [([], None), (["123", "456"], "123")])
def test_doi2pmid(id_list: list[str], expected: str | None) -> None:
    # PubMed returns an empty IdList for DOIs it doesn't index (FPBASE-73F)
    with (
        patch.object(entrez.Entrez, "esearch"),
        patch.object(
            entrez.Entrez, "read", return_value={"Count": str(len(id_list)), "IdList": id_list}
        ),
    ):
        assert entrez._doi2pmid("10.1134/s1068162017010034") == expected
