const { test, expect } = require('@playwright/test');

test('comparison renders missing, failed and verified evidence without invented passes', async ({ page }) => {
  await page.goto('/');
  await page.waitForFunction(() => typeof renderComparisonDashboard === 'function');
  const render = async (rule, ml, extra = {}) => page.evaluate(({ rule, ml, extra }) => {
    document.getElementById('resultsContainer').classList.remove('hidden');
    renderComparisonDashboard({ mode: 'dual_compare', comparison: {
      rule: { metrics: rule }, ml: { metrics: ml }, ...extra
    } });
  }, { rule, ml, extra });
  const row = label => page.locator('#comparisonMatrixBody tr').filter({ hasText: label });

  await render({}, { aa_identity: null, cai: null, gc_percent: null, type_iis_clean: null });
  for (const label of ['AA Translation Identity', 'GC Content', 'CAI Index', 'Type IIS Clearance']) {
    await expect(row(label).locator('td').nth(1)).toHaveText('Not evaluated');
    await expect(row(label).locator('td').nth(2)).toHaveText('Not evaluated');
    await expect(row(label).locator('td').nth(3)).toHaveText('Not evaluated');
  }
  await expect(page.locator('#comparisonMatrixBody')).not.toContainText('NaN');
  await expect(page.locator('#comparisonMatrixBody')).not.toContainText('Passed');

  await render({ aa_identity: 1, type_iis_clean: true, type_iis_site_count: 0, cai: 0, gc_percent: 0 },
    { aa_identity: 0.5, type_iis_clean: false, type_iis_site_count: 2, cai: 0.8, gc_percent: 40 });
  await expect(row('AA Translation Identity')).toContainText('100.00% Passed');
  await expect(row('AA Translation Identity')).toContainText('50.00% Failed');
  await expect(row('AA Translation Identity').locator('td').nth(3)).toHaveText('Failed');
  await expect(row('Type IIS Clearance')).toContainText('2 site(s) — Failed');
  await expect(row('CAI Index').locator('td').nth(1)).toHaveText('0.000');
  await expect(row('GC Content').locator('td').nth(1)).toHaveText('0.0%');

  await render({ aa_identity: 1, type_iis_clean: true }, { aa_identity: 1, type_iis_clean: true });
  await expect(row('AA Translation Identity').locator('td').nth(3)).toHaveText('Passed');
  await expect(row('Type IIS Clearance').locator('td').nth(3)).toHaveText('Passed');

  await render({ aa_identity: '1', type_iis_clean: 'true', cai: '0.9', gc_percent: -1 },
    { aa_identity: 1.1, type_iis_clean: true, type_iis_site_count: 2, cai: false, gc_percent: 101 });
  for (const label of ['AA Translation Identity', 'GC Content', 'CAI Index', 'Type IIS Clearance']) {
    await expect(row(label).locator('td').nth(3)).toHaveText('Not evaluated');
  }
  await render({}, { type_iis_clean: false });
  await expect(row('Type IIS Clearance').locator('td').nth(2)).toHaveText('Failed');
  await expect(row('Type IIS Clearance').locator('td').nth(3)).toHaveText('Failed');
});
