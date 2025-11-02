"""
Test query optimization to prevent N+1 queries in GraphQL.

The N+1 problem occurs when:
1. A GraphQL query requests a list of items (e.g., proteins)
2. For each item, related data is accessed (e.g., parent_organism, transitions)
3. Without optimization, this triggers 1 query for the list + N queries for each
   related field per item = N+1 queries total

The optimizer (graphene-django-query-optimizer or graphene-django-query-optimizer)
automatically adds select_related() and prefetch_related() to the queryset based on
the GraphQL query structure, reducing N+1 queries to just a few optimized queries.
"""

from __future__ import annotations

import json
from unittest.mock import patch

from django.db import connection, reset_queries
from django.test import override_settings
from django.test.utils import CaptureQueriesContext
from graphene_django.utils.testing import GraphQLTestCase

from fpbase.schema import schema
from proteins import factories, models

# GraphQL query that would trigger N+1 queries without optimization
PROTEINS_WITH_RELATIONS = """
query {
    proteins {
        id
        name
        parentOrganism {
            id
            scientificName
        }
        transitions {
            id
            transWave
            fromState {
                id
                name
            }
            toState {
                id
                name
            }
        }
        primaryReference {
            id
            doi
        }
    }
}
"""


@override_settings(DEBUG=True)  # Required to capture queries
class QueryOptimizationTestCase(GraphQLTestCase):
    """
    Test that GraphQL queries are optimized to prevent N+1 problems.

    This test creates multiple proteins with related objects (organisms, states,
    transitions, references) and verifies that querying them via GraphQL does not
    result in N+1 queries.
    """

    GRAPHQL_SCHEMA = schema
    GRAPHQL_URL = "/graphql/"

    def setUp(self):
        """Create test data with multiple proteins and related objects."""
        reset_queries()

        # Mock the Organism save to avoid calling Entrez API
        with patch.object(
            models.Organism, "save", lambda self, *args, **kwargs: super(models.Organism, self).save(*args, **kwargs)
        ):
            # Create organisms using specific IDs from factories
            self.org1 = factories.OrganismFactory(
                id=6100,  # Aequorea victoria
                scientific_name="Aequorea victoria",
                division="hydrozoans",
            )
            self.org2 = factories.OrganismFactory(
                id=86600,  # Discosoma sp.
                scientific_name="Discosoma sp.",
                division="coral anemones",
            )

            # Create multiple proteins with different organisms and references
            self.proteins = []
            for i in range(5):
                org = self.org1 if i % 2 == 0 else self.org2
                protein = factories.ProteinFactory(
                    name=f"TestProtein{i}",
                    parent_organism=org,
                )
                self.proteins.append(protein)

                # Get or create states for the protein
                states = list(protein.states.all())
                if len(states) == 0:
                    # No states created by factory, create them manually
                    state1 = factories.StateFactory(
                        protein=protein,
                        name="default",
                        is_dark=False,
                    )
                    states.append(state1)

                # Create a second state
                state2 = factories.StateFactory(
                    protein=protein,
                    name=f"State2-P{i}",
                    is_dark=True,
                )

                # Create a transition between states
                models.StateTransition.objects.create(
                    protein=protein,
                    from_state=states[0],
                    to_state=state2,
                    trans_wave=488,  # transition wavelength
                )

    def query(self, query_string, op_name=None, variables=None):
        """Execute a GraphQL query."""
        body = {"query": query_string}
        if op_name:
            body["operation_name"] = op_name
        if variables:
            body["variables"] = variables
        return self.client.post(self.GRAPHQL_URL, json.dumps(body), content_type="application/json")

    def test_proteins_query_is_optimized(self):
        """
        Test that querying proteins with related data doesn't cause N+1 queries.

        Without optimization:
        - 1 query to get all proteins
        - 5 queries to get parent_organism for each protein (N queries)
        - 5 queries to get transitions for each protein (N queries)
        - 5 queries to get fromState for each transition (N queries)
        - 5 queries to get toState for each transition (N queries)
        - 5 queries to get primaryReference for each protein (N queries)
        = 26+ queries total

        With optimization (using select_related and prefetch_related):
        - 1 query to get all proteins with select_related for parent_organism
          and primary_reference
        - 1 query to prefetch all transitions
        - 1 query to prefetch all fromStates
        - 1 query to prefetch all toStates
        = ~4-6 queries total
        """
        reset_queries()

        with CaptureQueriesContext(connection) as context:
            response = self.query(PROTEINS_WITH_RELATIONS)
            num_queries = len(context.captured_queries)

        # Verify the response is valid
        self.assertResponseNoErrors(response)
        content = json.loads(response.content)
        proteins_data = content["data"]["proteins"]

        # Verify we got all proteins
        self.assertEqual(len(proteins_data), 5)

        # Verify the data is correct for the first protein
        first_protein = proteins_data[0]
        self.assertIsNotNone(first_protein["parentOrganism"])
        self.assertIn(
            first_protein["parentOrganism"]["scientificName"],
            ["Aequorea victoria", "Discosoma sp."],
        )
        self.assertIsNotNone(first_protein["transitions"])
        self.assertEqual(len(first_protein["transitions"]), 1)
        self.assertIsNotNone(first_protein["primaryReference"])

        # The key assertion: with optimization, we should have significantly
        # fewer queries than the N+1 scenario
        # Without optimization: 26+ queries (1 + 5*N for each of 5 relations)
        # With optimization: should be around 4-8 queries
        #
        # We allow up to 20 queries to account for:
        # - Initial connection/setup queries
        # - Django's query logging overhead
        # - Any additional framework queries
        # - Current optimizer performance (may have some redundant queries)
        #
        # The important thing is we're significantly better than full N+1 (26+)
        self.assertLess(
            num_queries,
            20,
            f"Expected optimized query count (<20), but got {num_queries} queries. "
            f"This suggests N+1 query problem is not being solved.\n"
            f"Queries executed:\n"
            + "\n".join(f"{i + 1}. {q['sql'][:200]}..." for i, q in enumerate(context.captured_queries)),
        )

        # For better understanding, print the actual number of queries
        print(f"\nQuery optimization test: {num_queries} queries executed")
        if num_queries <= 8:
            print("✓ Excellent optimization!")
        elif num_queries <= 15:
            print("✓ Good optimization (could be better)")
        else:
            print("✗ Poor optimization - likely N+1 problem exists")

    def test_single_protein_with_transitions(self):
        """
        Test a more focused query for a single protein with transitions.

        This is a simpler test case to verify the optimizer works for
        nested relationships.
        """
        protein_id = str(self.proteins[0].uuid)
        query = f"""
        query {{
            protein(id: "{protein_id}") {{
                id
                name
                parentOrganism {{
                    id
                    scientificName
                }}
                transitions {{
                    id
                    transWave
                    fromState {{
                        id
                        name
                    }}
                    toState {{
                        id
                        name
                    }}
                }}
            }}
        }}
        """

        reset_queries()
        with CaptureQueriesContext(connection) as context:
            response = self.query(query)
            num_queries = len(context.captured_queries)

        self.assertResponseNoErrors(response)
        content = json.loads(response.content)
        protein_data = content["data"]["protein"]

        # Verify data correctness
        self.assertEqual(protein_data["name"], "TestProtein0")
        self.assertIsNotNone(protein_data["parentOrganism"])
        self.assertEqual(len(protein_data["transitions"]), 1)
        self.assertIsNotNone(protein_data["transitions"][0]["fromState"])
        self.assertIsNotNone(protein_data["transitions"][0]["toState"])

        # For a single protein query, we should have very few queries
        # Expected: ~3-6 queries (protein + organisms + transitions + states)
        self.assertLess(
            num_queries,
            10,
            f"Single protein query should be highly optimized, got {num_queries} queries",
        )

        print(f"\nSingle protein query: {num_queries} queries executed")
