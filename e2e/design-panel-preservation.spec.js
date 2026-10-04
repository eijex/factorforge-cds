const { test, expect } = require('@playwright/test');
const { execFileSync } = require('node:child_process');

const baselineJs = execFileSync('git', ['show', 'caa405c:web/js/app.js'], { encoding: 'utf8' });
const baselineHtml = execFileSync('git', ['show', 'caa405c:web/index.html'], { encoding: 'utf8' });

async function capture(context, baseline, preset, overrides) {
  const page = await context.newPage();
  let payload;
  if (baseline) {
    await page.route('**/js/app.js', route => route.fulfill({ contentType: 'text/javascript', body: baselineJs }));
    await page.route('http://127.0.0.1:8010/', route => route.fulfill({ contentType: 'text/html', body: baselineHtml }));
  }
  await page.route('**/api/optimize', async route => {
    if (route.request().method() === 'GET') {
      await route.fulfill({ json: { supported_objectives: ['dp_v2_1_1'], capabilities: {}, validation_checks: [] } });
    } else {
      payload = route.request().postDataJSON();
      await route.fulfill({ status: 400, json: { error: 'Synthetic contract capture only' } });
    }
  });
  await page.goto('/');
  await expect(page.locator('#sopProfileMeta')).not.toHaveText('');
  await page.locator('#sopPresetSelect').selectOption(preset);
  if (overrides) {
    await page.locator('#manualSopOverrides > summary').click();
    await page.locator('#optimizationSeed').fill('0');
    await page.locator('#optimizationSeed').dispatchEvent('change');
    await page.locator('#criterionCaiMode').selectOption('required');
    await page.locator('input[name="typeIisEnzyme"][value="SapI"]').check();
  }
  const inventory = await page.locator('input[id], select[id], textarea[id]').evaluateAll(nodes => nodes.map(n => n.id).sort());
  await page.locator('#sequenceInput').fill('MASTWQHLEKDPVNRSGYF');
  await expect(page.locator('#optimizeBtn')).toBeEnabled();
  await page.locator('#optimizeBtn').click();
  await expect.poll(() => payload).toBeTruthy();
  await page.close();
  return { payload, inventory };
}

for (const preset of ['preset_a', 'preset_b', 'preset_c', 'preset_d', 'preset_e']) {
  test(`existing ${preset} request and controls are preserved`, async ({ browser }) => {
    const before = await browser.newContext();
    const after = await browser.newContext();
    try {
      expect(await capture(after, false, preset, false)).toEqual(await capture(before, true, preset, false));
    } finally { await before.close(); await after.close(); }
  });
}

test('seed zero, review mode and restriction override requests are preserved', async ({ browser }) => {
  const before = await browser.newContext();
  const after = await browser.newContext();
  try {
    expect(await capture(after, false, 'preset_a', true)).toEqual(await capture(before, true, 'preset_a', true));
  } finally { await before.close(); await after.close(); }
});

for (const width of [390, 768, 1440]) {
  test(`disclosures retain accessible controls at ${width}px`, async ({ page }) => {
    await page.setViewportSize({ width, height: 900 });
    await page.goto('/');
    await expect(page.locator('#engineVersionSelect')).toBeHidden();
    await expect(page.locator('#toggleWatermark')).toBeHidden();
    await page.locator('#manualSopOverrides > summary').click();
    await expect(page.locator('#engineVersionSelect')).toBeVisible();
    await expect(page.locator('#useTemplate')).toBeVisible();
    await page.locator('#experimentalSettings > summary').click();
    await expect(page.locator('#toggleWatermark')).toBeVisible();
    const overflow = await page.evaluate(() => Array.from(document.querySelectorAll('main *')).filter(n => {
      const r = n.getBoundingClientRect();
      return r.width && r.right > innerWidth + 1;
    }).map(n => n.id || n.tagName).slice(0, 12));
    expect(overflow).toEqual([]);
  });
}
