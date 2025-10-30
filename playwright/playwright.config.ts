/// <reference types="node" />

import path from 'node:path';

import { defineConfig, devices } from '@playwright/test';

type ReporterTuple = [string] | [string, Record<string, unknown>];

const port = process.env.PLAYWRIGHT_DEV_SERVER_PORT ?? '8000';
const baseURL = process.env.PLAYWRIGHT_BASE_URL ?? `http://127.0.0.1:${port}`;
const webServerCommand =
  process.env.PLAYWRIGHT_WEB_SERVER_COMMAND ?? 'uv run backend/manage.py runserver 0.0.0.0:8000';
const skipWebServer = process.env.PLAYWRIGHT_SKIP_WEB_SERVER === '1';
const repoRoot = path.resolve(__dirname, '..');

const reporters: ReporterTuple[] = [['list']];
if (!process.env.CI) {
  reporters.push(['html', { open: 'never' }]);
}

export default defineConfig({
  testDir: './tests',
  timeout: 60_000,
  expect: {
    timeout: 10_000,
  },
  fullyParallel: true,
  forbidOnly: !!process.env.CI,
  retries: process.env.CI ? 2 : 0,
  workers: process.env.CI ? 2 : undefined,
  reporter: process.env.CI ? 'github' : reporters,
  use: {
    baseURL,
    trace: 'on-first-retry',
    screenshot: 'only-on-failure',
    video: 'retain-on-failure',
  },
  projects: [
    {
      name: 'chromium',
      use: { ...devices['Desktop Chrome'] },
    },
    {
      name: 'firefox',
      use: { ...devices['Desktop Firefox'] },
    },
    {
      name: 'webkit',
      use: { ...devices['Desktop Safari'] },
    },
  ],
  webServer: skipWebServer
    ? undefined
    : {
        command: webServerCommand,
        url: baseURL,
        reuseExistingServer: true,
        cwd: repoRoot,
        timeout: 120_000,
      },
});
