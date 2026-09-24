const { test, expect } = require('@playwright/test');

test('default SOP is applied and a custom upload persists locally', async ({ page }) => {
  await page.goto('/');
  await expect(page.locator('#sopProfileName')).toHaveText('Conservative Plant Expression Review Template');
  await expect(page.locator('#manualSopOverrides')).not.toHaveAttribute('open', '');
  await expect(page.locator('#sopUploadZone')).toContainText('.yaml, .yml, or .json');
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

test('YAML SOP upload is validated, applied, and persisted', async ({ page }) => {
  await page.goto('/');
  const yaml = `$schema: factorforge-sop-v1
profile_id: local_yaml_v1
sop_name: "Local YAML Review SOP"
version: 1.0.0
author: "Local user"
status: LOCAL
default_enforcement: IGNORE
unknown_rule_policy: ERROR
rules:
  rna.cryptic_splice.v1: WARNING
workflow:
  design_method: gc_target
  comparison_methods: [high_cai]
  sequence_requirements:
    reproducibility_seed: 42
    kozak_optimization: yes
    suppress_tpa: no
    type_iis_enzymes: [BsaI]
    custom_restriction_sites: []
    moclo_template: no
  review_policy:
    cai: PREFERRED
    overall_gc: REQUIRED
    local_gc: PREFERRED
    type_iis: REQUIRED
    repeats: PREFERRED
    homopolymers: PREFERRED
    forbidden_motifs: PREFERRED
`;
  await page.locator('#sopFileUpload').setInputFiles({ name: 'local-sop.yaml', mimeType: 'application/yaml', buffer: Buffer.from(yaml) });
  await expect(page.locator('#sopProfileName')).toHaveText('Local YAML Review SOP');
  await expect(page.locator('input[name="objective"][value="gc_target"]')).toBeChecked();
  await expect(page.locator('#toggleKozak')).toBeChecked();
  await page.reload();
  await expect(page.locator('#sopProfileName')).toHaveText('Local YAML Review SOP');
});

test('invalid SOP upload is rejected without replacing the active profile', async ({ page }) => {
  await page.goto('/');
  await page.locator('#sopFileUpload').setInputFiles({ name: 'bad.json', mimeType: 'application/json', buffer: Buffer.from('{"profile_id":"bad"}') });
  await expect(page.locator('#sopProfileName')).toHaveText('Conservative Plant Expression Review Template');
  await expect(page.locator('#toastContainer')).toContainText('SOP rejected');
});

test('default SOP summary is available without expanding manual controls', async ({ page }) => {
  await page.goto('/');
  await page.locator('#defaultSopSummary summary').click();
  await expect(page.locator('#defaultSopSummary')).toContainText('BsaI, BpiI, and BsmBI');
  await expect(page.locator('#manualSopOverrides')).not.toHaveAttribute('open', '');
});

test('active SOP downloads as human-readable YAML', async ({ page }) => {
  await page.goto('/');
  const downloadPromise = page.waitForEvent('download');
  await page.locator('#downloadSopTemplate').click();
  const download = await downloadPromise;
  expect(download.suggestedFilename()).toMatch(/\.yaml$/);
  const stream = await download.createReadStream();
  let text = '';
  for await (const chunk of stream) text += chunk.toString();
  expect(text).toContain('$schema: factorforge-sop-v1');
  expect(text).toContain('type_iis_enzymes: [BsaI, BpiI, BsmBI]');
});

test('bundled example SOP is directly downloadable and sharing stays explicit', async ({ page }) => {
  await page.goto('/');
  const example = page.locator('#downloadSopExample');
  await expect(example).toHaveAttribute('href', '/examples/factorforge-conservative-sop.yaml');
  await expect(example).toHaveAttribute('download', '');
  await expect(example).toHaveAttribute('title', /commented FactorForge conservative example/);
  await expect(page.locator('#downloadSopTemplate')).toHaveAttribute('title', /currently active/);
  await expect(page.locator('#uploadSopButton')).toHaveAttribute('title', /YAML or JSON/);
  await expect(page.getByRole('link', { name: /Share a public-safe SOP suggestion/ })).toHaveAttribute('href', /template=sop_template\.yml/);

  const response = await page.request.get('/examples/factorforge-conservative-sop.yaml');
  expect(response.ok()).toBeTruthy();
  const yaml = await response.text();
  expect(yaml).toContain('# This is an in-silico starting template');
  expect(yaml).toContain('$schema: factorforge-sop-v1');
  expect(yaml).toContain('unknown_rule_policy: ERROR');

  await page.locator('#sopFileUpload').setInputFiles('web/examples/factorforge-conservative-sop.yaml');
  await expect(page.locator('#sopProfileName')).toHaveText('FactorForge Conservative Plant Expression Review Template');
});
