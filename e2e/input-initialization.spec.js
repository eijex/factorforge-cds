const { test, expect } = require('@playwright/test');

test('input works while server metadata is still pending', async ({ page }) => {
  let release;
  const gate = new Promise(resolve => { release = resolve; });
  await page.route('**/api/optimize', async route => {
    await gate;
    await route.fulfill({ json: { capabilities: {}, validation_checks: [] } });
  });
  await page.goto('/', { waitUntil: 'domcontentloaded' });
  await page.locator('#sequenceInput').fill('MASTWQHLEKDPVNRSGYF');
  try {
    await expect(page.locator('#optimizeBtn')).toBeEnabled({ timeout: 2000 });
  } finally {
    release();
  }
});
