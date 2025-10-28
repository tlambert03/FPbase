/**
 * Apollo Client Configuration Validation Tests
 *
 * These tests verify that the Apollo Client configuration is correct,
 * particularly the possibleTypes configuration which is critical for
 * fragments on interface types to work.
 */

import { describe, expect, it } from 'vitest'
import introspectionQueryResultData from '../../fragmentTypes.json'
import { getPossibleTypes, validatePossibleTypes } from '../fixtures/apolloClient'

describe('Apollo Client possibleTypes Configuration', () => {
  it('has valid introspection data structure', () => {
    expect(introspectionQueryResultData).toBeDefined()
    expect(introspectionQueryResultData.__schema).toBeDefined()
    expect(introspectionQueryResultData.__schema.types).toBeInstanceOf(Array)
  })

  it('contains FluorophoreInterface in introspection data', () => {
    const fluorophoreType = introspectionQueryResultData.__schema.types.find(
      (type) => type.name === 'FluorophoreInterface'
    )

    expect(fluorophoreType).toBeDefined()
    expect(fluorophoreType.kind).toBe('INTERFACE')
    expect(fluorophoreType.possibleTypes).toBeInstanceOf(Array)
  })

  it('FluorophoreInterface has State and Dye as possible types', () => {
    const fluorophoreType = introspectionQueryResultData.__schema.types.find(
      (type) => type.name === 'FluorophoreInterface'
    )

    const typeNames = fluorophoreType.possibleTypes.map((t) => t.name)
    expect(typeNames).toContain('State')
    expect(typeNames).toContain('Dye')
  })

  it('contains SpectrumOwnerInterface in introspection data', () => {
    const ownerType = introspectionQueryResultData.__schema.types.find(
      (type) => type.name === 'SpectrumOwnerInterface'
    )

    expect(ownerType).toBeDefined()
    expect(ownerType.kind).toBe('INTERFACE')
    expect(ownerType.possibleTypes).toBeInstanceOf(Array)
  })

  it('SpectrumOwnerInterface has all expected possible types', () => {
    const ownerType = introspectionQueryResultData.__schema.types.find(
      (type) => type.name === 'SpectrumOwnerInterface'
    )

    const typeNames = ownerType.possibleTypes.map((t) => t.name)
    const expectedTypes = ['State', 'Dye', 'Filter', 'Light', 'Camera']

    expectedTypes.forEach((expectedType) => {
      expect(typeNames).toContain(expectedType)
    })
  })

  it('transforms introspection data to Apollo v3 possibleTypes format', () => {
    const possibleTypes = getPossibleTypes()

    // Should be an object, not the original introspection structure
    expect(possibleTypes).toBeTypeOf('object')
    expect(possibleTypes.__schema).toBeUndefined()

    // Should have interface names as keys
    expect(possibleTypes).toHaveProperty('FluorophoreInterface')
    expect(possibleTypes).toHaveProperty('SpectrumOwnerInterface')
  })

  it('FluorophoreInterface maps to array of type names', () => {
    const possibleTypes = getPossibleTypes()

    const fluorophoreTypes = possibleTypes.FluorophoreInterface
    expect(fluorophoreTypes).toBeInstanceOf(Array)
    expect(fluorophoreTypes).toContain('State')
    expect(fluorophoreTypes).toContain('Dye')

    // Should be strings, not objects
    fluorophoreTypes.forEach((type) => {
      expect(typeof type).toBe('string')
    })
  })

  it('SpectrumOwnerInterface maps to array of all owner type names', () => {
    const possibleTypes = getPossibleTypes()

    const ownerTypes = possibleTypes.SpectrumOwnerInterface
    expect(ownerTypes).toBeInstanceOf(Array)

    const expectedTypes = ['State', 'Dye', 'Filter', 'Light', 'Camera']
    expectedTypes.forEach((expectedType) => {
      expect(ownerTypes).toContain(expectedType)
    })
  })

  it('validatePossibleTypes passes for correct configuration', () => {
    // Should not throw
    expect(() => validatePossibleTypes()).not.toThrow()
    expect(validatePossibleTypes()).toBe(true)
  })

  it('getPossibleTypes matches production client transformation', () => {
    // This is the EXACT transformation used in production client.js
    const productionTransform = introspectionQueryResultData.__schema.types.reduce((acc, type) => {
      if (type.possibleTypes) {
        acc[type.name] = type.possibleTypes.map((t) => t.name)
      }
      return acc
    }, {})

    const testUtilityTransform = getPossibleTypes()

    // They must be identical
    expect(testUtilityTransform).toEqual(productionTransform)
  })
})

describe('Fragment Configuration Regression', () => {
  it('prevents regression: possibleTypes must be defined', () => {
    const possibleTypes = getPossibleTypes()

    // The bug was caused by possibleTypes being undefined
    // This test ensures that never happens again
    expect(possibleTypes).toBeDefined()
    expect(possibleTypes).not.toBeNull()
    expect(Object.keys(possibleTypes).length).toBeGreaterThan(0)
  })

  it('prevents regression: FluorophoreInterface must be in possibleTypes', () => {
    const possibleTypes = getPossibleTypes()

    // The bug manifested because fragments on FluorophoreInterface didn't work
    expect(possibleTypes).toHaveProperty('FluorophoreInterface')

    const fluorophoreTypes = possibleTypes.FluorophoreInterface
    expect(fluorophoreTypes).toBeInstanceOf(Array)
    expect(fluorophoreTypes.length).toBeGreaterThan(0)
  })

  it('prevents regression: State and Dye must implement FluorophoreInterface', () => {
    const possibleTypes = getPossibleTypes()

    const fluorophoreTypes = possibleTypes.FluorophoreInterface

    // If either State or Dye is missing, fragments won't apply to them
    expect(fluorophoreTypes).toContain('State')
    expect(fluorophoreTypes).toContain('Dye')
  })

  it('prevents regression: introspection format must be compatible', () => {
    // The bug occurred during Apollo 2 → 3 upgrade
    // This verifies the introspection structure is usable

    expect(introspectionQueryResultData.__schema).toBeDefined()
    expect(introspectionQueryResultData.__schema.types).toBeInstanceOf(Array)

    const hasInterfaces = introspectionQueryResultData.__schema.types.some(
      (type) => type.kind === 'INTERFACE' && type.possibleTypes
    )

    expect(hasInterfaces).toBe(true)
  })
})
