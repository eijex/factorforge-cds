const {test, expect} = require('@playwright/test');

test('unconnected comparisons stay off and disabled, without changing engine settings', async ({page}) => {
  await page.route('**/api/optimize', route => route.fulfill({json:{validation_checks:[],capabilities:{}}}));
  await page.goto('/');
  await expect(page.locator('#experimentalComparisonCard')).not.toHaveAttribute('open', '');
  await page.locator('#experimentalComparisonCard > summary').click();
  await expect(page.locator('#experimentalComparisonCard')).toContainText('not connected to web execution');
  for (const name of ['toggleArgSerDistribution','toggleSpliceScan','toggleFivePrimeBoundary']) {
    await expect(page.locator(`input[name="${name}"]`)).not.toBeChecked();
    await expect(page.locator(`input[name="${name}"]`)).toBeDisabled();
  }
  for (const name of ['argSerDistributionMode','spliceScanMode','fivePrimeBoundaryMode']) {
    await expect(page.locator(`select[name="${name}"]`)).toBeDisabled();
  }
  await expect(page.locator('select[name="spliceScanMode"]')).toHaveValue('warn_only');
  await expect(page.locator('select[name="fivePrimeBoundaryMode"]')).toHaveValue('free');
  await expect(page.locator('#experimentalComparisonCard')).not.toContainText('for translation efficiency');
  await page.reload();
  await expect(page.locator('input[name="toggleSpliceScan"]')).not.toBeChecked();
  await expect(page.locator('input[name="toggleFivePrimeBoundary"]')).not.toBeChecked();
});
