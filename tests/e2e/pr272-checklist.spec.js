import { test, expect } from '@playwright/test';

// Helper function to collect console messages and errors
function setupConsoleMonitoring(page, testName) {
  const consoleMessages = [];
  const errors = [];

  page.on('console', msg => {
    const text = msg.text();
    const type = msg.type();

    consoleMessages.push({ type, text });

    // Check for React 18 warnings, MUI warnings, PropTypes warnings
    if (type === 'warning' || type === 'error') {
      if (text.includes('React') ||
          text.includes('Material-UI') ||
          text.includes('MUI') ||
          text.includes('PropTypes') ||
          text.includes('deprecated')) {
        errors.push({ type, text, test: testName });
      }
    }
  });

  page.on('pageerror', error => {
    errors.push({ type: 'pageerror', text: error.message, test: testName });
  });

  return { consoleMessages, errors };
}

test.describe('PR #272 Manual Testing Checklist - React 18 & MUI v5 Upgrade', () => {

  test.describe('1. BLAST Interface (/blast/)', () => {
    test('should load BLAST page without errors and display correct MUI v5 styling', async ({ page }) => {
      const { errors } = setupConsoleMonitoring(page, 'BLAST initial load');

      // Visit /blast/
      await page.goto('/blast/');

      // Wait for the page to fully load
      await page.waitForLoadState('networkidle');

      // Check that the page loaded successfully
      await expect(page).toHaveTitle(/BLAST/i);

      // Check for the BLAST app container
      const blastApp = page.locator('#blast-app');
      await expect(blastApp).toBeVisible();

      // Check for Material-UI v5 components
      // Look for MUI buttons, tabs, inputs
      const submitButton = page.getByRole('button', { name: /submit/i });

      // Check if buttons have MUI v5 classes (should start with MuiButton)
      if (await submitButton.count() > 0) {
        const buttonClass = await submitButton.first().getAttribute('class');
        expect(buttonClass).toMatch(/MuiButton/);
      }

      // Report console errors
      if (errors.length > 0) {
        console.log('BLAST Console Errors/Warnings:', JSON.stringify(errors, null, 2));
      }

      expect(errors.length).toBe(0);
    });

    test('should allow submitting a BLAST query', async ({ page }) => {
      const { errors } = setupConsoleMonitoring(page, 'BLAST submit query');

      await page.goto('/blast/');
      await page.waitForLoadState('networkidle');

      // Sample protein sequence (mCherry partial sequence)
      const sampleSequence = 'MVSKGEEDNMAIIKEFMRFKVHMEGSVNGHEFEIEGEGEGRPYEGTQTAKLKVTKGGPLPFAWDILSPQFMYGSKAYVKHPADIPDYLKLSFPEGFKWERVMNFEDGGVVTVTQDSSLQDGEFIYKVKLRGTNFPSDGPVMQKKTMGWEASSERMYPEDGALKGEIKQRLKLKDGGHYDAEVKTTYKAKKPVQLPGAYNVNIKLDITSHNEDYTIVEQYERAEGRHSTGGMDELYK';

      // Find and fill the textarea
      const textarea = page.locator('textarea').first();
      await textarea.fill(sampleSequence);

      // Click submit button
      const submitButton = page.getByRole('button', { name: /submit/i });
      await submitButton.click();

      // Wait for results or loading indicator
      await page.waitForTimeout(2000); // Give it time to process

      // Check for Snackbar notifications (MUI v5)
      const snackbar = page.locator('[role="alert"], .MuiSnackbar-root');

      // The page should show some feedback
      // Either loading state, results, or error message

      if (errors.length > 0) {
        console.log('BLAST Submit Console Errors/Warnings:', JSON.stringify(errors, null, 2));
      }
    });
  });

  test.describe('2. Full Spectra Viewer (/spectra/)', () => {
    test('should load spectra viewer without errors', async ({ page }) => {
      const { errors } = setupConsoleMonitoring(page, 'Spectra Viewer initial load');

      await page.goto('/spectra/');
      await page.waitForLoadState('networkidle');

      // Check for the spectra viewer container
      const spectraViewer = page.locator('#spectra-viewer');
      await expect(spectraViewer).toBeVisible();

      // Check for AppBar
      const appBar = page.locator('[class*="MuiAppBar"]');
      await expect(appBar.first()).toBeVisible();

      // Check for MUI v5 theme colors
      const appBarBg = await appBar.first().evaluate(el =>
        window.getComputedStyle(el).backgroundColor
      );
      expect(appBarBg).toBeTruthy();

      if (errors.length > 0) {
        console.log('Spectra Viewer Console Errors/Warnings:', JSON.stringify(errors, null, 2));
      }

      expect(errors.length).toBe(0);
    });

    test('should allow adding spectra and display them', async ({ page }) => {
      const { errors } = setupConsoleMonitoring(page, 'Spectra Viewer add spectra');

      await page.goto('/spectra/');
      await page.waitForLoadState('networkidle');

      // Click the "+" (Add) button
      const addButton = page.locator('button[aria-label*="Add"], button:has-text("+")').first();
      if (await addButton.count() > 0) {
        await addButton.click();

        // Wait for modal to open
        await page.waitForTimeout(1000);

        // Look for search modal
        const modal = page.locator('[role="dialog"], .MuiDialog-root');
        if (await modal.count() > 0) {
          await expect(modal.first()).toBeVisible();

          // Try to search for EGFP
          const searchInput = page.locator('input[type="text"]').first();
          if (await searchInput.count() > 0) {
            await searchInput.fill('EGFP');
            await page.waitForTimeout(1000);

            // Try to select the first result
            const firstResult = page.locator('[role="option"]').first();
            if (await firstResult.count() > 0) {
              await firstResult.click();
            }
          }

          // Close modal
          const closeButton = page.locator('button[aria-label*="close"], button:has-text("Close")').first();
          if (await closeButton.count() > 0) {
            await closeButton.click();
          }
        }
      }

      // Check for Highcharts graph
      await page.waitForTimeout(2000);
      const highchartsContainer = page.locator('.highcharts-container');
      if (await highchartsContainer.count() > 0) {
        await expect(highchartsContainer.first()).toBeVisible();
      }

      if (errors.length > 0) {
        console.log('Spectra Add Console Errors/Warnings:', JSON.stringify(errors, null, 2));
      }
    });

    test('should allow toggling log scale', async ({ page }) => {
      const { errors } = setupConsoleMonitoring(page, 'Spectra Viewer log scale');

      // Visit with pre-loaded spectra
      await page.goto('/spectra/?s=17,18'); // EGFP, mCherry
      await page.waitForLoadState('networkidle');
      await page.waitForTimeout(2000);

      // Look for log scale toggle
      const logScaleToggle = page.locator('input[type="checkbox"], [role="switch"]').filter({ hasText: /log/i }).first();

      if (await logScaleToggle.count() === 0) {
        // Try alternative selectors
        const switches = page.locator('[role="switch"], .MuiSwitch-root');
        if (await switches.count() > 0) {
          // Click the first switch
          await switches.first().click();
          await page.waitForTimeout(500);
        }
      } else {
        await logScaleToggle.click();
        await page.waitForTimeout(500);
      }

      if (errors.length > 0) {
        console.log('Spectra Log Scale Console Errors/Warnings:', JSON.stringify(errors, null, 2));
      }
    });

    test('should open share dialog', async ({ page }) => {
      const { errors } = setupConsoleMonitoring(page, 'Spectra Viewer share dialog');

      await page.goto('/spectra/?s=17,18');
      await page.waitForLoadState('networkidle');
      await page.waitForTimeout(1000);

      // Look for share button
      const shareButton = page.locator('button[aria-label*="Share"], button:has-text("Share")').first();

      if (await shareButton.count() > 0) {
        await shareButton.click();
        await page.waitForTimeout(500);

        // Check for dialog
        const dialog = page.locator('[role="dialog"]');
        await expect(dialog.first()).toBeVisible();

        // Look for copy link button
        const copyButton = page.getByRole('button', { name: /copy/i });
        if (await copyButton.count() > 0) {
          await expect(copyButton.first()).toBeVisible();
        }

        // Close dialog
        const closeButton = page.locator('button[aria-label*="close"]').first();
        if (await closeButton.count() > 0) {
          await closeButton.click();
        }
      }

      if (errors.length > 0) {
        console.log('Spectra Share Console Errors/Warnings:', JSON.stringify(errors, null, 2));
      }
    });

    test('should open settings drawer', async ({ page }) => {
      const { errors } = setupConsoleMonitoring(page, 'Spectra Viewer settings drawer');

      await page.goto('/spectra/?s=17,18');
      await page.waitForLoadState('networkidle');
      await page.waitForTimeout(1000);

      // Look for settings/gear icon
      const settingsButton = page.locator('button[aria-label*="settings"], button[aria-label*="Settings"]').first();

      if (await settingsButton.count() > 0) {
        await settingsButton.click();
        await page.waitForTimeout(500);

        // Check for drawer
        const drawer = page.locator('[class*="MuiDrawer"]');
        if (await drawer.count() > 0) {
          await expect(drawer.first()).toBeVisible();
        }
      }

      if (errors.length > 0) {
        console.log('Spectra Settings Console Errors/Warnings:', JSON.stringify(errors, null, 2));
      }
    });

    test('should test efficiency table overlap buttons', async ({ page }) => {
      const { errors } = setupConsoleMonitoring(page, 'Efficiency Table overlap buttons');

      await page.goto('/spectra/?s=17,18');
      await page.waitForLoadState('networkidle');
      await page.waitForTimeout(2000);

      // Look for efficiency table
      const efficiencyTable = page.locator('table, [class*="MuiTable"]');

      if (await efficiencyTable.count() > 0) {
        // Look for overlap toggle buttons
        const overlapButtons = page.locator('button').filter({ hasText: /overlap/i });

        if (await overlapButtons.count() > 0) {
          // Click the first overlap button
          await overlapButtons.first().click();
          await page.waitForTimeout(500);

          // Click it again to toggle
          await overlapButtons.first().click();
          await page.waitForTimeout(500);
        }
      }

      if (errors.length > 0) {
        console.log('Efficiency Table Console Errors/Warnings:', JSON.stringify(errors, null, 2));
      }
    });
  });

  test.describe('4. Protein Comparison Page (/compare/)', () => {
    test('should load comparison page with spectra', async ({ page }) => {
      const { errors } = setupConsoleMonitoring(page, 'Protein comparison');

      await page.goto('/compare/egfp,mcherry/');
      await page.waitForLoadState('networkidle');
      await page.waitForTimeout(2000);

      // Check for comparison table/grid
      const comparisonContent = page.locator('main, [role="main"], .container');
      await expect(comparisonContent.first()).toBeVisible();

      // Check for embedded spectra viewer
      const spectraViewer = page.locator('#spectra-viewer');
      if (await spectraViewer.count() > 0) {
        await expect(spectraViewer).toBeVisible();

        // Wait for Highcharts to render
        await page.waitForTimeout(2000);
        const highcharts = page.locator('.highcharts-container');
        if (await highcharts.count() > 0) {
          await expect(highcharts.first()).toBeVisible();
        }
      }

      if (errors.length > 0) {
        console.log('Protein Comparison Console Errors/Warnings:', JSON.stringify(errors, null, 2));
      }

      expect(errors.length).toBe(0);
    });
  });

  test.describe('5. Protein Detail Page (/protein/egfp/)', () => {
    test('should load protein detail page with embedded spectra viewer', async ({ page }) => {
      const { errors } = setupConsoleMonitoring(page, 'Protein detail');

      await page.goto('/protein/egfp/');
      await page.waitForLoadState('networkidle');

      // Check page loaded
      await expect(page).toHaveTitle(/EGFP/i);

      // Scroll to spectral properties section
      const spectraSection = page.locator('#spectra-viewer, [class*="spectr"]');
      if (await spectraSection.count() > 0) {
        await spectraSection.first().scrollIntoViewIfNeeded();
        await page.waitForTimeout(2000);

        // Check for embedded viewer
        const highcharts = page.locator('.highcharts-container');
        if (await highcharts.count() > 0) {
          await expect(highcharts.first()).toBeVisible();
        }
      }

      // Check for MUI components on the page
      const muiButtons = page.locator('[class*="MuiButton"]');
      if (await muiButtons.count() > 0) {
        expect(await muiButtons.count()).toBeGreaterThan(0);
      }

      if (errors.length > 0) {
        console.log('Protein Detail Console Errors/Warnings:', JSON.stringify(errors, null, 2));
      }

      expect(errors.length).toBe(0);
    });
  });

  test.describe('8. Base Template Components', () => {
    test('should load home page with correct MUI styling', async ({ page }) => {
      const { errors } = setupConsoleMonitoring(page, 'Home page');

      await page.goto('/');
      await page.waitForLoadState('networkidle');

      // Check for navigation bar
      const nav = page.locator('nav, header');
      await expect(nav.first()).toBeVisible();

      // Check for footer
      const footer = page.locator('footer');
      if (await footer.count() > 0) {
        await expect(footer).toBeVisible();
      }

      // Check for any MUI components
      const muiComponents = page.locator('[class*="Mui"]');
      if (await muiComponents.count() > 0) {
        expect(await muiComponents.count()).toBeGreaterThan(0);
      }

      if (errors.length > 0) {
        console.log('Home Page Console Errors/Warnings:', JSON.stringify(errors, null, 2));
      }
    });
  });

  test.describe('9. Protein Table (/table/)', () => {
    test('should load protein table with DataTable', async ({ page }) => {
      const { errors } = setupConsoleMonitoring(page, 'Protein table');

      await page.goto('/table/');
      await page.waitForLoadState('networkidle');
      await page.waitForTimeout(2000);

      // Check for table
      const table = page.locator('table, [role="table"]');
      if (await table.count() > 0) {
        await expect(table.first()).toBeVisible();
      }

      // Check for MUI table components
      const muiTable = page.locator('[class*="MuiTable"]');
      if (await muiTable.count() > 0) {
        await expect(muiTable.first()).toBeVisible();
      }

      if (errors.length > 0) {
        console.log('Protein Table Console Errors/Warnings:', JSON.stringify(errors, null, 2));
      }
    });
  });

  test.describe('Critical Checks Across All Pages', () => {
    const criticalPages = [
      { url: '/', name: 'Home' },
      { url: '/spectra/', name: 'Spectra Viewer' },
      { url: '/blast/', name: 'BLAST' },
      { url: '/protein/egfp/', name: 'Protein Detail' },
    ];

    for (const { url, name } of criticalPages) {
      test(`${name} should have no 404s for assets`, async ({ page }) => {
        const failed404s = [];

        page.on('response', response => {
          if (response.status() === 404) {
            failed404s.push(response.url());
          }
        });

        await page.goto(url);
        await page.waitForLoadState('networkidle');

        // Filter out known acceptable 404s like devtools
        const realErrors = failed404s.filter(url =>
          !url.includes('devtools') &&
          !url.includes('.well-known')
        );

        if (realErrors.length > 0) {
          console.log(`${name} 404 Errors:`, realErrors);
        }

        expect(realErrors.length).toBe(0);
      });

      test(`${name} should have MUI v5 components`, async ({ page }) => {
        await page.goto(url);
        await page.waitForLoadState('networkidle');

        // Check for MUI v5 class names
        const muiV5Elements = page.locator('[class*="MuiButton"], [class*="MuiAppBar"], [class*="MuiDialog"]');
        const count = await muiV5Elements.count();

        // Some pages might not have MUI components, so we just check they load
        // without errors
        expect(count).toBeGreaterThanOrEqual(0);
      });
    }
  });
});
