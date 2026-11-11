const { test, expect } = require('@playwright/test');

test('public 1DP/3DP/4DP availability and custom SOP AGPL-3.0 terms remain consistent', async ({ page }, testInfo) => {
  await page.route('**/api/optimize', route => route.fulfill({
    json: { capabilities: { direct_structure: { available: true, status: 'experimental' }, local_4dp: { available: true, fixtures: false } } }
  }));
  await page.goto('/');
  await expect(page.locator('input[value="baseline_1dp"]')).toBeChecked();
  await expect(page.locator('input[value="paired_3dp"]')).toBeEnabled();
  await expect(page.locator('input[value="paired_4dp"]')).toBeEnabled();
  await expect(page.locator('#policyScopeNotice')).toContainText('not an approved laboratory SOP');
  await expect(page.locator('#policyScopeNotice')).toContainText('GNU AGPL-3.0');
  await expect(page.locator('body')).toContainText('eijex.lab@gmail.com');
  await expect(page.locator('body')).toContainText('shared with Eijex and FactorForge');
  await page.screenshot({ path: testInfo.outputPath('public-policy.png'), fullPage: true });
});

test('explicit API pause disables 3DP and 4DP cleanly', async ({ page }) => {
  await page.route('**/api/optimize', route => route.fulfill({
    json: { capabilities: { direct_structure: { available: false, status: 'paused' }, local_4dp: { available: false, fixtures: false } } }
  }));
  await page.goto('/');
  await expect(page.locator('input[value="baseline_1dp"]')).toBeChecked();
  await expect(page.locator('input[value="paired_3dp"]')).toBeDisabled();
  await expect(page.locator('input[value="paired_4dp"]')).toBeDisabled();
  await page.evaluate(() => {
    const input = document.querySelector('input[value="paired_4dp"]');
    input.disabled = false;
    input.checked = true;
    input.dispatchEvent(new Event('change'));
  });
  await expect(page.locator('input[value="baseline_1dp"]')).toBeChecked();
  await expect(page.locator('#structure4dpInputs')).toBeHidden();
});
