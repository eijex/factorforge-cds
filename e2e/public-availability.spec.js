const { test, expect } = require('@playwright/test');

test('public research hold and data-use terms remain consistent', async ({ page }, testInfo) => {
  await page.route('**/api/optimize', route => route.fulfill({
    json: { capabilities: { direct_structure: { available: false, status: 'paused' } } }
  }));
  await page.goto('/');
  await expect(page.locator('input[value="baseline_1dp"]')).toBeChecked();
  await expect(page.locator('input[value="paired_3dp"]')).toBeDisabled();
  await expect(page.locator('#policyScopeNotice')).toContainText('not a laboratory-approved SOP');
  await expect(page.locator('body')).toContainText('does not automatically transfer ownership');
  await expect(page.locator('body')).toContainText('requires prior explicit consent');
  await page.screenshot({ path: testInfo.outputPath('public-policy.png'), fullPage: true });
  await page.evaluate(() => {
    const input = document.querySelector('input[value="paired_3dp"]');
    input.disabled = false;
    input.checked = true;
    input.dispatchEvent(new Event('change'));
  });
  await expect(page.locator('input[value="baseline_1dp"]')).toBeChecked();
  await expect(page.locator('#structureWebInputs')).toBeHidden();
});
