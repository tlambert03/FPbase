import bs4
import pytest
from django import forms

from proteins.forms.spectrum import SpectrumForm
from proteins.models.state import State
from tomcomplete import AutocompleteSelect


@pytest.mark.django_db
def test_select_widget():
    field = forms.ModelChoiceField(
        label="Protein",
        queryset=State.objects.select_related("protein"),
        widget=AutocompleteSelect(
            url="proteins:state-autocomplete",
            attrs={"id": "id_test_select"},
        ),
    )
    NAME = "test_select"
    html = field.widget.render(name=NAME, value="default")
    select = bs4.BeautifulSoup(html, "html.parser").find("select")
    assert select is not None
    assert select["name"] == NAME
    assert select.has_attr("id")

    first_option = select.find("option")
    assert first_option is not None
    assert first_option.text == "---------"  # django default for empty choice


@pytest.mark.django_db
def test_spectrum_form_autocomplete():
    form = SpectrumForm()
    print(str(form))
