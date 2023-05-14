import json

from graphene_django.utils.testing import GraphQLTestCase

from proteins import models

from ..schema import schema

OPTICAL_CONFIG = """
query OpticalConfig($id: Int!) {
    opticalConfig(id: $id) {
        id
        name
        microscope {
            id
            name
        }
        filters {
            id
            path
            reflects
            spectrum {
            id
            }
        }
        light {
            id
            spectrum {
            id
            }
        }
        camera {
            id
            spectrum {
            id
            }
        }
        laser
        comments
    }
}
"""

SPECTRUM = """
query Spectrum($id: Int!) {
    spectrum(id: $id) {
        id
        data
        category
        color
        subtype
        owner {
            slug
            name
            id
            ... on State {
                ...FluorophoreParts
            }
            ... on Dye {
                ...FluorophoreParts
            }
        }
    }
}

fragment FluorophoreParts on FluorophoreInterface {
    qy
    extCoeff
    twopPeakgm
    exMax
    emMax
}
"""


class SpectraQueriesTestCase(GraphQLTestCase):
    GRAPHQL_SCHEMA = schema
    GRAPHQL_URL = "/graphql/"  # default is '/graphql' ... which gives a redirect

    def setUp(self):
        self.microscope = models.Microscope.objects.create()
        self.protein = models.Protein.objects.create(name="test")
        self.optical_config = models.OpticalConfig.objects.create(microscope=self.microscope)
        self.dye = models.Dye.objects.get_or_create(name="test-dye")[0]
        self.spectrum = models.Spectrum.objects.get_or_create(
            category=models.Spectrum.DYE,
            subtype=models.Spectrum.EM,
            owner_dye=self.dye,
            data=[[0, 1], [1, 1]],
        )[0]

    def query(self, query, op_name=None, input_data=None, variables=None):
        body = {"query": query}
        if op_name:
            body["operation_name"] = op_name
        if variables:
            body["variables"] = variables
        if input_data:
            body["variables"] = {"input": input_data}
            if variables in body and isinstance(body, dict):
                body["variables"]["input"] = input_data
            else:
                body["variables"] = {"input": input_data}
        return self.client.post(self.GRAPHQL_URL, json.dumps(body), content_type="application/json")

    def test_optical_configs(self):
        response = self.query(
            """
            query   {
                opticalConfigs {
                    id
                    name
                    comments
                    microscope {
                        id
                        name
                    }
                }
            }
            """,
            op_name="opticalConfig",
        )
        # content = json.loads(response.content)

        self.assertResponseNoErrors(response)

    def test_optical_config_works(self):
        response = self.query(OPTICAL_CONFIG, variables={"id": self.optical_config.id})
        content = json.loads(response.content)
        self.assertResponseNoErrors(response)
        self.assertEqual(content["data"]["opticalConfig"]["microscope"]["id"], self.microscope.id)

    def test_optical_config_breaks(self):
        response = self.query(OPTICAL_CONFIG, variables={"id": 99999})
        self.assertResponseHasErrors(response)

    def test_spectrum(self):
        response = self.query(SPECTRUM, op_name="Spectrum", variables={"id": self.spectrum.id})
        self.assertResponseNoErrors(response)
        content = json.loads(response.content)
        self.assertEqual(content["data"]["spectrum"]["owner"]["id"], str(self.dye.id))
        self.assertEqual(
            content["data"]["spectrum"]["category"].upper(),
            str(self.spectrum.category.upper()),
        )

    def test_spectra(self):
        response = self.query(
            """
            {
                spectra {
                    id
                    category
                    subtype
                    owner {
                        name
                        slug
                        url
                    }
                }
            }
            """
        )
        self.assertResponseNoErrors(response)
        content = json.loads(response.content)
        last_spectrum = content["data"]["spectra"][-1]
        self.assertEqual(last_spectrum["id"], str(self.spectrum.id))
        self.assertEqual(last_spectrum["owner"]["name"], self.dye.name)

    def test_protein(self):
        response = self.query(
            """
            query Protein($id: String!) {
                protein(id: $id){
                    name
                    id
                    slug
                    aliases
                    seqValidated
                    seq
                    seqComment
                    pdb
                    genbank
                    uniprot
                    ipgId
                    mw
                    references {
                        doi
                    }
                }
            }
            """,
            variables={"id": self.protein.uuid},
        )
        self.assertResponseNoErrors(response)
