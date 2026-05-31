const { test, expect } = require('@playwright/test');

const SAMPLE_PROTEIN = 'MSKGEELFTGVVPILVELD';
const MOCK_DNA = 'ATGTCCAAGGGCGAGGAGCTGTTCACCGGCGTGGTGCCCATCCTGGTGGAGCTGGAC';

async function openApp(page) {
  await page.addInitScript(() => {
    window.Chart = class {
      destroy() {}
    };
  });
  await page.goto('/');
}

test('loads the main web UI', async ({ page }) => {
  await openApp(page);

  await expect(page.locator('#sequenceInput')).toBeVisible();
  await expect(page.getByRole('heading', { name: '⚙️ Optimization Settings' })).toBeVisible();
  await expect(page.locator('#optimizeBtn')).toBeVisible();
});

test('updates sequence metadata for protein input', async ({ page }) => {
  await openApp(page);

  await page.locator('#sequenceInput').fill(SAMPLE_PROTEIN);

  await expect(page.locator('#inputTypeBadge')).toHaveText('Protein');
  await expect(page.locator('#inputLenBadge')).toHaveText(`${SAMPLE_PROTEIN.length} aa`);
  await expect(page.locator('#sequencePreview')).toContainText(SAMPLE_PROTEIN);
});

test('BY-2 host disables feasibility_best and selects a profile fallback', async ({ page }) => {
  await openApp(page);

  await page.getByText('Tobacco BY-2', { exact: true }).click();

  await expect(page.locator('input[name="host"][value="by2"]')).toBeChecked();
  await expect(page.locator('input[name="objective"][value="feasibility_best"]')).toBeDisabled();
  await expect(page.locator('input[name="objective"][value="high_cai"]')).toBeChecked();
});

test('optimization payload includes host and renders host_profile', async ({ page }) => {
  let requestBody;
  await page.route('**/api/optimize', async route => {
    requestBody = route.request().postDataJSON();
    await route.fulfill({
      status: 200,
      contentType: 'application/json',
      body: JSON.stringify({
        success: true,
        optimized_sequence: MOCK_DNA,
        original_length: SAMPLE_PROTEIN.length,
        optimized_length: MOCK_DNA.length,
        metrics: {
          cai: 0.91,
          gc_percent: 58.3,
          polya_signals: 0,
          length: MOCK_DNA.length
        },
        profile: 'high_cai',
        host_profile: 'by2',
        validation: {
          input_type: 'protein',
          polya: 'PASS',
          moclo: 'UNCHECKED',
          gc: 'PASS'
        }
      })
    });
  });
  await openApp(page);

  await page.locator('#sequenceInput').fill(SAMPLE_PROTEIN);
  await page.getByText('Tobacco BY-2', { exact: true }).click();
  await page.locator('#optimizeBtn').click();

  await expect.poll(() => requestBody).toMatchObject({
    sequence: SAMPLE_PROTEIN,
    host: 'by2',
    profile: 'high_cai'
  });
  await expect(page.locator('#hostProfileValue')).toContainText('by2');
  await expect(page.locator('#optimizedSequence')).toContainText(MOCK_DNA.slice(0, 20));
});

test('clear input resets preview and sequence badges', async ({ page }) => {
  await openApp(page);

  await page.locator('#sequenceInput').fill(SAMPLE_PROTEIN);
  await expect(page.locator('#sequencePreview')).toBeVisible();

  await page.locator('#clearBtn').click();

  await expect(page.locator('#sequenceInput')).toHaveValue('');
  await expect(page.locator('#previewContainer')).toBeHidden();
  await expect(page.locator('#inputLenBadge')).toHaveText('0 bp');
});
