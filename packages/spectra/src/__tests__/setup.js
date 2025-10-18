import '@testing-library/jest-dom'
import { expect, afterEach } from 'vitest'
import { cleanup } from '@testing-library/react'

// Cleanup after each test
afterEach(() => {
  cleanup()
})

// Mock Highcharts globally to avoid initialization issues in tests
global.Highcharts = {
  charts: [],
  addEvent: () => {},
  removeEvent: () => {},
}
