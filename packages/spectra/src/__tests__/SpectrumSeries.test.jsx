import { render } from "@testing-library/react"
import Highcharts from "highcharts"
import { List } from "immutable"
import { HighchartsChart, HighchartsProvider, XAxis, YAxis } from "react-jsx-highcharts"
import { describe, expect, it, vi } from "vitest"
import SpectrumSeries from "../Components/SpectraViewer/SpectrumSeries"
import { mockSpectrumEGFP_EM, mockSpectrumEGFP_EX } from "./mockData"

// Mock Apollo Client
const mockApolloClient = {
  query: vi.fn(),
}

vi.mock("@apollo/client", () => ({
  useApolloClient: () => mockApolloClient,
}))

// Wrapper component that provides the necessary Highcharts context
const TestWrapper = ({ children }) => {
  return (
    <HighchartsProvider Highcharts={Highcharts}>
      <HighchartsChart>
        <XAxis id="xAxis">
          <YAxis id="yAxis">{children}</YAxis>
        </XAxis>
      </HighchartsChart>
    </HighchartsProvider>
  )
}

describe("SpectrumSeries", () => {
  it("renders without crashing with excitation spectrum", () => {
    const { container } = render(
      <TestWrapper>
        <SpectrumSeries
          spectrum={mockSpectrumEGFP_EX}
          areaFill={true}
          palette="wavelength"
          ownerIndex={0}
          ownerInfo={{}}
          visible={true}
        />
      </TestWrapper>
    )
    expect(container).toBeTruthy()
  })

  it("renders without crashing with emission spectrum", () => {
    const { container } = render(
      <TestWrapper>
        <SpectrumSeries
          spectrum={mockSpectrumEGFP_EM}
          areaFill={true}
          palette="wavelength"
          ownerIndex={0}
          ownerInfo={{}}
          visible={true}
        />
      </TestWrapper>
    )
    expect(container).toBeTruthy()
  })

  it("handles Immutable.js List data correctly", () => {
    // Create spectrum with Immutable.js List data (simulates real app behavior)
    const spectrumWithImmutableData = {
      ...mockSpectrumEGFP_EX,
      data: List(mockSpectrumEGFP_EX.data.map((point) => List(point))),
    }

    const { container } = render(
      <TestWrapper>
        <SpectrumSeries
          spectrum={spectrumWithImmutableData}
          areaFill={true}
          palette="wavelength"
          ownerIndex={0}
          ownerInfo={{}}
          visible={true}
        />
      </TestWrapper>
    )

    // Should render without errors
    expect(container).toBeTruthy()
  })

  it("renders with scaleEC enabled for excitation spectrum", () => {
    const { container } = render(
      <TestWrapper>
        <SpectrumSeries
          spectrum={mockSpectrumEGFP_EX}
          areaFill={true}
          scaleEC={true}
          palette="wavelength"
          ownerIndex={0}
          ownerInfo={{}}
          visible={true}
        />
      </TestWrapper>
    )
    expect(container).toBeTruthy()
  })

  it("renders with scaleQY enabled for emission spectrum", () => {
    const { container } = render(
      <TestWrapper>
        <SpectrumSeries
          spectrum={mockSpectrumEGFP_EM}
          areaFill={true}
          scaleQY={true}
          palette="wavelength"
          ownerIndex={0}
          ownerInfo={{}}
          visible={true}
        />
      </TestWrapper>
    )
    expect(container).toBeTruthy()
  })

  it("renders with logScale enabled", () => {
    const { container } = render(
      <TestWrapper>
        <SpectrumSeries
          spectrum={mockSpectrumEGFP_EX}
          areaFill={true}
          logScale={true}
          palette="wavelength"
          ownerIndex={0}
          ownerInfo={{}}
          visible={true}
        />
      </TestWrapper>
    )
    expect(container).toBeTruthy()
  })

  it("renders with inverted data", () => {
    const { container } = render(
      <TestWrapper>
        <SpectrumSeries
          spectrum={mockSpectrumEGFP_EX}
          areaFill={true}
          inverted={true}
          palette="wavelength"
          ownerIndex={0}
          ownerInfo={{}}
          visible={true}
        />
      </TestWrapper>
    )
    expect(container).toBeTruthy()
  })

  it("renders multiple series together", () => {
    const { container } = render(
      <TestWrapper>
        <SpectrumSeries
          spectrum={mockSpectrumEGFP_EX}
          areaFill={true}
          palette="wavelength"
          ownerIndex={0}
          ownerInfo={{}}
          visible={true}
        />
        <SpectrumSeries
          spectrum={mockSpectrumEGFP_EM}
          areaFill={true}
          palette="wavelength"
          ownerIndex={0}
          ownerInfo={{}}
          visible={true}
        />
      </TestWrapper>
    )
    expect(container).toBeTruthy()
  })
})
