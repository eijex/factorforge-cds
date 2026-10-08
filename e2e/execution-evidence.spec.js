const { test, expect } = require('@playwright/test');

test('execution trace does not invent passed checks or engine provenance', async ({ page }) => {
  await page.route('**/api/optimize', route => route.fulfill({ json: { supported_objectives: [], validation_checks: [] } }));
  await page.goto('/');
  await page.evaluate(() => renderExecutionEvidence({}));
  const grid = page.locator('#traceStagesGrid');
  await expect(grid).not.toContainText('Complete');
  await expect(grid).not.toContainText('screens passed');
  await expect(grid).toContainText('Not reported');
  await expect(grid).toContainText('Not verified here');
  await page.evaluate(() => renderExecutionEvidence({ design_contract: { engine_id: 'synthetic-engine', engine_version: '0.0.1' } }));
  await expect(grid).toContainText('synthetic-engine · 0.0.1');
  await expect(grid).not.toContainText('DP v2.1.1');
});
