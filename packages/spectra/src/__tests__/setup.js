import '@testing-library/jest-dom'
import { expect, beforeAll, afterAll, afterEach } from 'vitest'
import { cleanup } from '@testing-library/react'
import { setupServer } from 'msw/node'
import { handlers } from './msw/handlers'

// Set up MSW server for API mocking
const server = setupServer(...handlers)

// Start server before all tests
beforeAll(() => {
  server.listen({
    onUnhandledRequest: 'warn' // Warn on unhandled requests instead of error
  })

  // Suppress Highcharts warnings in test environment (no WebGL support in happy-dom)
  const originalConsoleError = console.error
  console.error = (...args) => {
    // Suppress Highcharts warning #26 (Boost module WebGL fallback)
    if (typeof args[0] === 'string' && args[0].includes('Highcharts warning')) {
      return
    }
    originalConsoleError.apply(console, args)
  }
})

// Reset handlers after each test
afterEach(() => {
  server.resetHandlers()
  cleanup()
})

// Close server after all tests
afterAll(() => {
  server.close()
})

// Export server for tests that need to add custom handlers
export { server }

// Mock Highcharts globally to avoid initialization issues in tests
global.Highcharts = {
  charts: [],
  addEvent: () => {},
  removeEvent: () => {},
}
