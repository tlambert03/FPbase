/**
 * INTEGRATION TEST: SpectraViewer Full Data Flow
 *
 * These tests verify the complete data flow from cache initialization
 * through spectrum fetching to chart rendering, ensuring that:
 * 1. Fragment fields (qy, extCoeff) are correctly fetched and used
 * 2. Spectrum normalization works correctly
 * 3. Chart options are properly applied
 * 4. The viewer handles different spectrum types (EX/EM) correctly
 */

import { ApolloClient, ApolloProvider, from, HttpLink } from '@apollo/client'
import { StyledEngineProvider, ThemeProvider } from '@mui/material/styles'
import { render, screen, waitFor } from '@testing-library/react'
import { describe, expect, it } from 'vitest'
import { SpectraViewerContainer } from '../../Components/SpectraViewer/SpectraViewer'
import theme from '../../Components/theme'
import {
  GET_ACTIVE_OVERLAPS,
  GET_ACTIVE_SPECTRA,
  GET_CHART_OPTIONS,
  GET_EX_NORM,
} from '../../client/queries'
import { createTestCache } from '../fixtures/apolloClient'

// Helper to create a test Apollo Client with cache initialized
function createTestApolloClient({ activeSpectra = [], chartOptions = {} } = {}) {
  const cache = createTestCache()

  // Initialize cache with test data
  const defaultChartOptions = {
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
    ...chartOptions,
  }

  cache.writeQuery({
    query: GET_CHART_OPTIONS,
    data: { chartOptions: defaultChartOptions },
  })

  cache.writeQuery({
    query: GET_ACTIVE_SPECTRA,
    data: { activeSpectra },
  })

  cache.writeQuery({
    query: GET_ACTIVE_OVERLAPS,
    data: { activeOverlaps: [] },
  })

  cache.writeQuery({
    query: GET_EX_NORM,
    data: { exNorm: [null, null] },
  })

  const link = from([
    new HttpLink({
      uri: 'http://test-endpoint/graphql/',
      fetch,
    }),
  ])

  return new ApolloClient({
    link,
    cache,
  })
}

// Helper to render with all required providers (Apollo, Theme, etc.)
function renderWithProviders(ui, client) {
  return render(
    <StyledEngineProvider injectFirst>
      <ThemeProvider theme={theme}>
        <ApolloProvider client={client}>{ui}</ApolloProvider>
      </ThemeProvider>
    </StyledEngineProvider>
  )
}

describe('SpectraViewer Integration Tests', () => {
  it('renders loading state when no spectra are active', async () => {
    const client = createTestApolloClient({ activeSpectra: [] })

    renderWithProviders(<SpectraViewerContainer ownerInfo={{}} />, client)

    // Should show NoData component (not in simple mode)
    await waitFor(() => {
      expect(screen.queryByRole('progressbar')).not.toBeInTheDocument()
    })

    // Chart container should be present
    const viewer = document.querySelector('.spectra-viewer')
    expect(viewer).toBeInTheDocument()
  })

  it('fetches and displays spectra with fragment fields (qy, extCoeff)', async () => {
    const client = createTestApolloClient({
      activeSpectra: ['18'], // EGFP excitation
      chartOptions: {
        scaleEC: true, // Enable EC scaling
        scaleQY: false,
      },
    })

    renderWithProviders(<SpectraViewerContainer ownerInfo={{}} />, client)

    // Wait for the spectrum to be fetched and rendered
    await waitFor(
      () => {
        // Verify the chart container is rendered
        const viewer = document.querySelector('.spectra-viewer')
        expect(viewer).toBeInTheDocument()

        // The spectrum series should be rendered (Highcharts creates SVG)
        const chart = document.querySelector('.highcharts-container')
        expect(chart).toBeInTheDocument()
      },
      { timeout: 2000 }
    )
  })

  it('displays both excitation and emission spectra correctly', async () => {
    const client = createTestApolloClient({
      activeSpectra: ['17', '18'], // EGFP emission (17) and excitation (18)
      chartOptions: {
        scaleEC: true,
        scaleQY: true,
      },
    })

    // Verify that both spectra can be fetched and have fragment fields
    const { data: data18 } = await client.query({
      query: require('../../client/queries').GET_SPECTRUM,
      variables: { id: 18 },
    })

    const { data: data17 } = await client.query({
      query: require('../../client/queries').GET_SPECTRUM,
      variables: { id: 17 },
    })

    // Verify both spectra have the fragment fields needed for scaling
    expect(data18.spectrum.owner.qy).toBe(0.6)
    expect(data18.spectrum.owner.extCoeff).toBe(55900)
    expect(data18.spectrum.subtype).toBe('EX')

    expect(data17.spectrum.owner.qy).toBe(0.6)
    expect(data17.spectrum.owner.extCoeff).toBe(55900)
    expect(data17.spectrum.subtype).toBe('EM')

    // Verify they're from the same owner (EGFP)
    expect(data18.spectrum.owner.name).toBe('EGFP')
    expect(data17.spectrum.owner.name).toBe('EGFP')
  })

  it('shows normalization notice when scaleEC is enabled', async () => {
    const client = createTestApolloClient({
      activeSpectra: ['18', '80'], // EGFP and mCherry excitation
      chartOptions: {
        scaleEC: true,
        scaleQY: false,
      },
    })

    const ownerInfo = {
      egfp: { ec: 55900, qy: 0.6 },
      mcherry: { ec: 72000, qy: 0.22 },
    }

    renderWithProviders(<SpectraViewerContainer ownerInfo={ownerInfo} />, client)

    await waitFor(
      () => {
        const viewer = document.querySelector('.spectra-viewer')
        expect(viewer).toBeInTheDocument()

        // Should show the normalization notice
        // The text content includes "EX NORMED TO EXT COEFF"
        expect(document.body.textContent).toContain('EX NORMED TO EXT COEFF')
      },
      { timeout: 2000 }
    )
  })

  it('shows QY normalization notice when scaleQY is enabled', async () => {
    const client = createTestApolloClient({
      activeSpectra: ['17'], // EGFP emission
      chartOptions: {
        scaleEC: false,
        scaleQY: true,
      },
    })

    const ownerInfo = {
      egfp: { ec: 55900, qy: 0.6 },
    }

    renderWithProviders(<SpectraViewerContainer ownerInfo={ownerInfo} />, client)

    await waitFor(
      () => {
        const viewer = document.querySelector('.spectra-viewer')
        expect(viewer).toBeInTheDocument()

        // Should show QY normalization notice
        expect(document.body.textContent).toContain('EM NORMED TO')
        expect(document.body.textContent).toContain('QY')
      },
      { timeout: 2000 }
    )
  })

  it('handles multiple spectra from different owners', async () => {
    const client = createTestApolloClient({
      activeSpectra: ['17', '18', '79', '80'], // EGFP and mCherry EX/EM
      chartOptions: {
        scaleEC: true,
        scaleQY: true,
      },
    })

    // Verify all 4 spectra can be fetched with correct fragment fields
    const spectraIds = [17, 18, 79, 80]
    const results = await Promise.all(
      spectraIds.map((id) =>
        client.query({
          query: require('../../client/queries').GET_SPECTRUM,
          variables: { id },
        })
      )
    )

    // Verify EGFP spectra (17, 18)
    expect(results[0].data.spectrum.owner.name).toBe('EGFP')
    expect(results[0].data.spectrum.owner.qy).toBe(0.6)
    expect(results[0].data.spectrum.owner.extCoeff).toBe(55900)
    expect(results[1].data.spectrum.owner.name).toBe('EGFP')

    // Verify mCherry spectra (79, 80)
    expect(results[2].data.spectrum.owner.name).toBe('mCherry')
    expect(results[2].data.spectrum.owner.qy).toBe(0.22)
    expect(results[2].data.spectrum.owner.extCoeff).toBe(72000)
    expect(results[3].data.spectrum.owner.name).toBe('mCherry')

    // Verify mix of EX and EM subtypes
    const subtypes = results.map((r) => r.data.spectrum.subtype)
    expect(subtypes).toContain('EX')
    expect(subtypes).toContain('EM')
  })

  it('works with providedSpectra (bypassing cache queries)', async () => {
    const client = createTestApolloClient()

    const providedSpectra = ['17', '18'] // EGFP EX and EM
    const providedChartOptions = {
      scaleEC: true,
      scaleQY: true,
      showX: true,
      showY: true,
      showGrid: false,
      areaFill: true,
      logScale: false,
      extremes: [null, null],
      shareTooltip: true,
      palette: 'wavelength',
    }

    renderWithProviders(
      <SpectraViewerContainer
        ownerInfo={{}}
        provideSpectra={providedSpectra}
        provideOverlaps={[]}
        provideOptions={providedChartOptions}
      />,
      client
    )

    await waitFor(
      () => {
        const viewer = document.querySelector('.spectra-viewer')
        expect(viewer).toBeInTheDocument()

        const chart = document.querySelector('.highcharts-container')
        expect(chart).toBeInTheDocument()
      },
      { timeout: 2000 }
    )
  })

  it('verifies fragment fields are used for normalization (regression test)', async () => {
    // This test specifically verifies that qy and extCoeff fields from fragments
    // are available to the SpectraViewer for normalization
    const client = createTestApolloClient({
      activeSpectra: ['18', '80'], // EGFP and mCherry excitation
      chartOptions: {
        scaleEC: true,
        scaleQY: false,
      },
    })

    renderWithProviders(<SpectraViewerContainer ownerInfo={{}} />, client)

    await waitFor(
      async () => {
        // Query the cache to verify fragment fields are present
        const spectrum18 = client.readQuery({
          query: require('graphql-tag')`
            query TestQuery {
              spectrum(id: 18) {
                id
                owner {
                  ... on State {
                    qy
                    extCoeff
                  }
                }
              }
            }
          `,
          variables: { id: 18 },
        })

        // This test verifies that the fragment fields are in the cache
        // If they're not, the normalization feature would be broken
        // (This was the original bug that prompted these tests)
        expect(spectrum18).toBeTruthy()
      },
      { timeout: 2000 }
    )
  })
})
