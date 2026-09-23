const { test, expect } = require('@playwright/test');

test('default SOP is applied and a custom upload persists locally', async ({ page }) => {
  await page.goto('/');
  await expect(page.locator('#sopProfileName')).toHaveText('Conservative Plant Expression Review Template');
  await expect(page.locator('input[name="typeIisEnzyme"][value="BsaI"]')).toBeChecked();
  await expect(page.locator('input[name="typeIisEnzyme"][value="BpiI"]')).toBeChecked();
  await expect(page.locator('#criterionTypeIisMode')).toHaveValue('required');

  const custom = {
    $schema: 'factorforge-sop-v1', profile_id: 'local_lab_v1', sop_name: 'Local Lab Review SOP',
    version: '1.0.0', author: 'Local user', status: 'LOCAL', default_enforcement: 'IGNORE',
    unknown_rule_policy: 'ERROR', rules: { 'rna.cryptic_splice.v1': 'WARNING' },
    workflow: {
      design_method: 'gc_target', comparison_methods: ['high_cai'],
      sequence_requirements: { reproducibility_seed: 42, kozak_optimization: true, suppress_tpa: false, type_iis_enzymes: ['BsaI'], custom_restriction_sites: [], moclo_template: false },
      review_policy: { cai: 'PREFERRED', overall_gc: 'REQUIRED', local_gc: 'PREFERRED', type_iis: 'REQUIRED', repeats: 'PREFERRED', homopolymers: 'PREFERRED', forbidden_motifs: 'PREFERRED' }
    }
  };
  await page.locator('#sopFileUpload').setInputFiles({ name: 'local-sop.json', mimeType: 'application/json', buffer: Buffer.from(JSON.stringify(custom)) });
  await expect(page.locator('#sopProfileName')).toHaveText('Local Lab Review SOP');
  await expect(page.locator('input[name="objective"][value="gc_target"]')).toBeChecked();
  await expect(page.locator('#toggleKozak')).toBeChecked();
  await page.reload();
  await expect(page.locator('#sopProfileName')).toHaveText('Local Lab Review SOP');
  await page.locator('#resetSopProfile').click();
  await expect(page.locator('#sopProfileName')).toHaveText('Conservative Plant Expression Review Template');
});

test('invalid SOP upload is rejected without replacing the active profile', async ({ page }) => {
  await page.goto('/');
  await page.locator('#sopFileUpload').setInputFiles({ name: 'bad.json', mimeType: 'application/json', buffer: Buffer.from('{"profile_id":"bad"}') });
  await expect(page.locator('#sopProfileName')).toHaveText('Conservative Plant Expression Review Template');
  await expect(page.locator('#toastContainer')).toContainText('SOP rejected');
});
