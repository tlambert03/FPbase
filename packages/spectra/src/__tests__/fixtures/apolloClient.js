/**
 * Test utilities for creating Apollo Client instances
 * Ensures tests use the same configuration as production
 */

import { InMemoryCache } from '@apollo/client'
import introspectionQueryResultData from '../../fragmentTypes.json'

/**
 * Transform Apollo v2 introspection format to Apollo v3 possibleTypes
 * This is the SAME transformation used in production client.js
 */
export function getPossibleTypes() {
  return introspectionQueryResultData.__schema.types.reduce((acc, type) => {
    if (type.possibleTypes) {
      acc[type.name] = type.possibleTypes.map((t) => t.name)
    }
    return acc
  }, {})
}

/**
 * Create a test cache with the same configuration as production
 * This ensures tests catch possibleTypes configuration bugs
 */
export function createTestCache() {
  const possibleTypes = getPossibleTypes()

  return new InMemoryCache({
    possibleTypes,
    // Note: canonizeResults option removed - deprecated in Apollo 3.14+
    typePolicies: {
      Query: {
        fields: {
          chartOptions: {
            merge(existing, incoming) {
              return { ...existing, ...incoming }
            },
          },
        },
      },
    },
  })
}

/**
 * Verify possibleTypes configuration is valid
 * This can be used in tests to ensure the configuration is correct
 */
export function validatePossibleTypes() {
  const possibleTypes = getPossibleTypes()

  // Critical: FluorophoreInterface must map to State and Dye
  const fluorophore = possibleTypes.FluorophoreInterface
  if (!fluorophore) {
    throw new Error('FluorophoreInterface not found in possibleTypes')
  }
  if (!fluorophore.includes('State') || !fluorophore.includes('Dye')) {
    throw new Error('FluorophoreInterface must include State and Dye types')
  }

  // Critical: SpectrumOwnerInterface must include all owner types
  const owner = possibleTypes.SpectrumOwnerInterface
  if (!owner) {
    throw new Error('SpectrumOwnerInterface not found in possibleTypes')
  }
  const expectedTypes = ['State', 'Dye', 'Filter', 'Light', 'Camera']
  for (const type of expectedTypes) {
    if (!owner.includes(type)) {
      throw new Error(`SpectrumOwnerInterface must include ${type} type`)
    }
  }

  return true
}
