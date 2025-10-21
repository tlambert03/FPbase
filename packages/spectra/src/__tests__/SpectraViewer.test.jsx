import { describe, it, expect } from 'vitest'
import { BaseSpectraViewer, XAxisWithRange, SpectraViewerContainer } from '../Components/SpectraViewer/SpectraViewer'

describe('SpectraViewer exports', () => {
  it('exports BaseSpectraViewer component', () => {
    expect(BaseSpectraViewer).toBeDefined()
    expect(BaseSpectraViewer).toBeTruthy()
  })

  it('exports XAxisWithRange component', () => {
    expect(XAxisWithRange).toBeDefined()
    expect(XAxisWithRange).toBeTruthy()
  })

  it('exports SpectraViewerContainer component', () => {
    expect(SpectraViewerContainer).toBeDefined()
    expect(SpectraViewerContainer).toBeTruthy()
  })
})
