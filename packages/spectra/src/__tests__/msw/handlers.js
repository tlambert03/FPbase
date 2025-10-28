/**
 * MSW (Mock Service Worker) handlers for GraphQL queries
 * These intercept network requests during tests and return mock data
 */

import { graphql, HttpResponse } from 'msw'
import { getSpectrumById } from '../fixtures/spectra'

// Create GraphQL handler for FPbase GraphQL endpoint
const fpbaseGraphQL = graphql.link('http://test-endpoint/graphql/')

export const handlers = [
  // GET_SPECTRUM query handler
  fpbaseGraphQL.query('Spectrum', ({ variables }) => {
    const { id } = variables
    const spectrum = getSpectrumById(Number(id))

    if (!spectrum) {
      return HttpResponse.json({
        errors: [{ message: `Spectrum with id ${id} not found` }],
      })
    }

    return HttpResponse.json({
      data: {
        spectrum,
      },
    })
  }),

  // Handle batch spectrum queries (if needed)
  fpbaseGraphQL.query('BatchSpectra', ({ variables }) => {
    const results = {}

    // BatchSpectra query uses dynamic field names like spectrum_17, spectrum_18
    Object.keys(variables).forEach((key) => {
      const match = key.match(/^spectrum_(\d+)$/)
      if (match) {
        const id = Number(match[1])
        const spectrum = getSpectrumById(id)
        if (spectrum) {
          results[key] = spectrum
        }
      }
    })

    return HttpResponse.json({
      data: results,
    })
  }),

  // GET_ACTIVE_SPECTRA query handler (client-side only, but can be mocked for integration tests)
  fpbaseGraphQL.query('ActiveSpectra', () => {
    return HttpResponse.json({
      data: {
        activeSpectra: ['17', '18', '79', '80'],
      },
    })
  }),

  // GET_CHART_OPTIONS query handler
  fpbaseGraphQL.query('ChartOptions', () => {
    return HttpResponse.json({
      data: {
        chartOptions: {
          __typename: 'ChartOptions',
          showY: true,
          showX: true,
          showGrid: false,
          areaFill: true,
          logScale: false,
          scaleEC: true,
          scaleQY: true,
          extremes: [null, null],
          shareTooltip: true,
          palette: 'wavelength',
        },
      },
    })
  }),
]

/**
 * Create a custom handler for specific test scenarios
 * Useful for testing error states or specific data combinations
 */
export function createSpectrumHandler(id, customData) {
  return fpbaseGraphQL.query('Spectrum', ({ variables }) => {
    if (variables.id === id) {
      return HttpResponse.json({
        data: {
          spectrum: customData,
        },
      })
    }
    // Fallback to default behavior
    const spectrum = getSpectrumById(Number(variables.id))
    return HttpResponse.json({
      data: {
        spectrum,
      },
    })
  })
}

/**
 * Create an error handler for testing error states
 */
export function createErrorHandler(errorMessage = 'GraphQL error') {
  return fpbaseGraphQL.query('Spectrum', () => {
    return HttpResponse.json({
      errors: [{ message: errorMessage }],
    })
  })
}
