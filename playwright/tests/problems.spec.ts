import { test, expect } from '@playwright/test';

const escapeRegex = (value: string): string => value.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');

const pages = [
  { path: '/problems/', heading: /Fixing Problems in FPbase/i },
  { path: '/problems/gaps/', heading: /Data gaps in FPbase/i },
  { path: '/problems/inconsistencies/', heading: /Possible database inconsistencies/i },
];

test.describe('Problem Pages', () => {
  for (const { path, heading } of pages) {
    test(`renders ${path}`, async ({ page }) => {
      await page.goto(path);
      await expect(page).toHaveURL(new RegExp(`${escapeRegex(path)}$`));
      await expect(page.getByRole('heading', { name: heading })).toBeVisible();
    });
  }
});
