import { test, expect } from '@playwright/test';

test.describe('Homepage', () => {
  test('shows navigation and search', async ({ page }) => {
    await page.goto('/');

    await expect(page.getByRole('navigation')).toBeVisible();
    await expect(page.getByRole('link', { name: 'FPbase' })).toBeVisible();

    const searchInput = page.getByPlaceholder('Search');
    await expect(searchInput).toBeVisible();
    await expect(searchInput).toBeEnabled();

    const mainContent = page.locator('main, #content, .main-content').first();
    await expect(mainContent).toBeVisible();
  });
});
