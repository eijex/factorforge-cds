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

test('discloses codon reference policy as current default plus non-selector note', async ({ page }) => {
  await openApp(page);

  const policy = page.locator('#codonReferencePolicy');
  await expect(policy).toBeVisible();
  await expect(policy).toContainText('Current default: NbeV1.1 HC CDS-derived');
  await expect(policy).toContainText('Selected');

  const packagedAssets = page.locator('#packagedReferenceAssets');
  await expect(packagedAssets).toContainText('not shown as public product choices');
  await expect(packagedAssets).toContainText('provenance, reproducibility, and controlled internal sensitivity analysis');
  await expect(packagedAssets).not.toContainText('Legacy Kazusa/SGN composite');
  await expect(packagedAssets).not.toContainText('NbeV1.1 all-CDS');
  await expect(packagedAssets).not.toContainText('QLD183 v103 CDS-derived');
  await expect(packagedAssets).not.toContainText('Tobacco BY-2 packaged table');
  await expect(policy.locator('input, select, button')).toHaveCount(0);
});

test('keeps non-default design objectives collapsed until requested', async ({ page }) => {
  await openApp(page);

  const objectives = page.locator('#designObjectivePolicy');
  await expect(objectives).toContainText('Feasibility Best');
  await expect(page.locator('input[name="objective"][value="feasibility_best"]')).toBeChecked();

  const implemented = page.locator('#implementedObjectives');
  const experimental = page.locator('#experimentalObjectives');
  await expect(implemented).not.toHaveAttribute('open', '');
  await expect(experimental).not.toHaveAttribute('open', '');
  await expect(implemented.getByText('High CAI')).toBeHidden();
  await expect(experimental.getByText("5' Ramp")).toBeHidden();

  await implemented.locator('summary').click();
  await expect(implemented).toHaveAttribute('open', '');
  await expect(implemented).toContainText('High CAI');
  await expect(implemented).toContainText('GC Target');
  await expect(implemented).toContainText('Assembly Friendly');

  await experimental.locator('summary').click();
  await expect(experimental).toHaveAttribute('open', '');
  await expect(experimental).toContainText("5' Ramp");
  await expect(experimental).toContainText('Viral Delivery');
  await expect(page.locator('input[name="objective"][value="ramp"]')).toBeDisabled();
  await expect(page.locator('input[name="objective"][value="viral_delivery"]')).toBeDisabled();
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
  await expect(page.locator('input[name="objective"][value="high_cai"]')).toBeDisabled();
  await expect(page.locator('input[name="objective"][value="gc_target"]')).toBeChecked();
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
        profile: 'gc_target',
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
    profile: 'gc_target'
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

async function fillAndOptimize(page) {
  await page.route('**/api/optimize', async route => {
    await route.fulfill({
      status: 200,
      contentType: 'application/json',
      body: JSON.stringify({
        success: true,
        optimized_sequence: MOCK_DNA,
        original_length: SAMPLE_PROTEIN.length,
        optimized_length: MOCK_DNA.length,
        metrics: { cai: 0.91, gc_percent: 58.3, polya_signals: 0, length: MOCK_DNA.length },
        profile: 'gc_target',
        host_profile: 'nbenthamiana',
        validation: { input_type: 'protein', polya: 'PASS', moclo: 'UNCHECKED', gc: 'PASS' }
      })
    });
  });
  await openApp(page);
  await page.locator('#sequenceInput').fill(SAMPLE_PROTEIN);
  await page.locator('#optimizeBtn').click();
  await expect(page.locator('#esmatlasFoldLink')).toHaveAttribute('href', /esmatlas\.com/);
}

test('third-party structure linkout is gated behind a consent modal (Cancel)', async ({ page }) => {
  await fillAndOptimize(page);

  await expect(page.locator('#linkoutConsentModal')).toHaveClass(/hidden/);

  await page.locator('#esmatlasFoldLink').click();
  await expect(page.locator('#linkoutConsentModal')).not.toHaveClass(/hidden/);
  await expect(page.locator('#linkoutConsentBody')).toContainText('Meta Platforms, Inc.');
  await expect(page.locator('#linkoutConsentBody')).not.toContainText('EvolutionaryScale');

  await page.locator('#linkoutConsentCancel').click();
  await expect(page.locator('#linkoutConsentModal')).toHaveClass(/hidden/);
});

test('third-party structure linkout consent (Continue) keeps the original AlphaFold DB URL', async ({ page }) => {
  await fillAndOptimize(page);

  const expectedHref = await page.locator('#alphafoldLink').getAttribute('href');
  await page.locator('#alphafoldLink').click();
  await expect(page.locator('#linkoutConsentModal')).not.toHaveClass(/hidden/);
  await expect(page.locator('#linkoutConsentBody')).toContainText('EMBL-EBI');
  await expect(page.locator('#linkoutConsentContinue')).toHaveAttribute('href', expectedHref);

  await page.locator('#linkoutConsentCancel').click();
});

test('structure linkout buttons are no-ops before any optimization result', async ({ page }) => {
  await openApp(page);

  // The buttons live inside #resultsContainer, which stays hidden until a
  // result renders — so this state is normally unreachable by a real click.
  // dispatchEvent('click') exercises the JS guard directly
  // (getAttribute('href') === '#') so a future markup change can't silently
  // remove that protection, without requiring the element to be visible.
  await expect(page.locator('#alphafoldLink')).toHaveAttribute('href', '#');
  await page.locator('#alphafoldLink').dispatchEvent('click');
  await expect(page.locator('#linkoutConsentModal')).toHaveClass(/hidden/);

  await expect(page.locator('#esmatlasFoldLink')).toHaveAttribute('href', '#');
  await page.locator('#esmatlasFoldLink').dispatchEvent('click');
  await expect(page.locator('#linkoutConsentModal')).toHaveClass(/hidden/);
});
