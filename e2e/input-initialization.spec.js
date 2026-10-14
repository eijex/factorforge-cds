const { test, expect } = require('@playwright/test');

// Retained custom SOP code is tested behind its disabled production feature flag.
test.beforeEach(async ({page}) => {
    const source = require('node:fs').readFileSync(require('node:path').join(__dirname, '../web/js/app.js'), 'utf8');
    await page.route('**/js/app.js', route => route.fulfill({contentType:'text/javascript', body:source.replace('const ENABLE_CUSTOM_SOP = false;', 'const ENABLE_CUSTOM_SOP = true;')}));
});

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
