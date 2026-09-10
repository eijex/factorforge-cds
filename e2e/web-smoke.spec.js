const { test, expect } = require('@playwright/test');

const SAMPLE_PROTEIN = 'MSKGEELFTGVVPILVELD';
const MOCK_DNA = 'ATGTCCAAGGGCGAGGAGCTGTTCACCGGCGTGGTGCCCATCCTGGTGGAGCTGGAC';
const REVIEW_ROWS = [
  { criterion: 'cai', mode: 'preferred', observed: 0.91, threshold: 0.8, result: 'PASS' },
  { criterion: 'overall_gc', mode: 'preferred', observed: 45, threshold: '40-47%', result: 'PASS' },
  { criterion: 'type_iis', mode: 'required', observed: 0, threshold: 0, result: 'PASS' },
];

function reviewResponse(overrides = {}) {
  return {
    success: true,
    optimized_sequence: MOCK_DNA,
    original_length: SAMPLE_PROTEIN.length,
    optimized_length: MOCK_DNA.length,
    input_type: 'protein',
    metrics: {
      cai: 0.91,
      gc_percent: 45,
      polya_signals: 0,
      length: MOCK_DNA.length,
      mfe_kcal_mol: null,
      mfe_status: 'not_computed',
      mfe_status_reason: 'missing_dependency',
      mfe_used: false,
      requested_gc_min_percent: 40,
      requested_gc_max_percent: 47,
    },
    profile: 'feasibility_best',
    host_profile: 'nbenthamiana',
    automated_decision: 'PASS',
    decision_summary: { required_failure_count: 0, preferred_warning_count: 0, explanation: 'All active acceptance criteria passed.' },
    qc_decision_matrix: REVIEW_ROWS,
    acceptance_criteria_snapshot: { cai: { mode: 'preferred', minimum: 0.8 } },
    validation: { input_type: 'protein', polya: 'PASS', moclo: 'PASS', gc: 'PASS' },
    constraint_report: { aa_identity: 1 },
    construct_id: 'CF-TEST-272',
    result_identifier: 'ff-result-272',
    created_at: '2026-09-09T00:00:00Z',
    product_version: '3.4.6',
    codon_reference_id: 'NbeV1.1-HC',
    reference_policy_version: '1.0',
    gc_reference_band: '40-47%',
    provenance: {
      input_sequence_hash: 'sha256:input-272',
      output_cds_hash: 'sha256:output-272',
      parameter_hash: 'sha256:params-272',
    },
    cds_design: { engine: 'factorforge_cds', objective: 'feasibility_best', product_version: '3.4.6' },
    ...overrides,
  };
}

async function mockOptimization(page, response) {
  await page.route('**/api/optimize', async route => {
    if (route.request().method() === 'GET') {
      await route.fulfill({ status: 200, contentType: 'application/json', body: JSON.stringify({ capabilities: {}, validation_checks: [] }) });
      return;
    }
    await route.fulfill({ status: 200, contentType: 'application/json', body: JSON.stringify(response) });
  });
}

async function openApp(page) {
  const pageErrors = [];
  page.on('pageerror', error => pageErrors.push(error.message));
  await page.addInitScript(() => {
    window.Chart = class {
      destroy() {}
    };
  });
  await page.goto('/');
  await page.waitForLoadState('networkidle');
  expect(pageErrors).toEqual([]);
}

test('loads the main web UI', async ({ page }) => {
  await openApp(page);

  await expect(page.locator('#sequenceInput')).toBeVisible();
  await expect(page.getByRole('heading', { name: '⚙️ Design Brief' })).toBeVisible();
  await expect(page.locator('#optimizeBtn')).toBeVisible();
  await expect(page.locator('#resultsPanel')).toBeVisible();
  await expect(page.locator('#emptyState')).toContainText('Awaiting a sequence');
  await expect(page.locator('#engineSelector')).toBeHidden();
  await expect(page.locator('#appliedPolicySummary')).toContainText('recommended feasibility design');
});

test('opens release notes and toggles dark mode', async ({ page }) => {
  await openApp(page);

  await page.locator('#themeToggle').click();
  await expect(page.locator('html')).toHaveClass(/dark/);

  await page.locator('#changelogBtn').click();
  await expect(page.locator('#changelogModal')).toBeVisible();
  await expect(page.locator('#changelogModal')).toContainText('v3.4.5');
  await page.locator('#closeModal').click();
  await expect(page.locator('#changelogModal')).toBeHidden();
});

test('keeps immutable codon-reference provenance out of the primary design form', async ({ page }) => {
  await openApp(page);

  const policy = page.locator('#codonReferencePolicy');
  await expect(policy).toBeHidden();
  await expect(policy).toContainText('Current default: NbeV1.1 HC CDS-derived');

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
  await expect(experimental).toBeHidden();
  await expect(implemented.getByText('High CAI')).toBeHidden();
  await expect(experimental.getByText("5' Ramp")).toBeHidden();

  await implemented.locator('summary').click();
  await expect(implemented).toHaveAttribute('open', '');
  await expect(implemented).toContainText('High CAI');
  await expect(implemented).toContainText('GC Target');
  await expect(implemented).toContainText('Assembly Friendly');

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

test('shows CDS design review controls and rejects multi-FASTA input', async ({ page }) => {
  await openApp(page);

  const acceptanceCriteria = page.locator('#acceptanceCriteria');
  await expect(acceptanceCriteria).toBeVisible();
  await expect(acceptanceCriteria).not.toHaveAttribute('open', '');
  await expect(page.locator('#criterionCaiMode')).toBeHidden();
  await acceptanceCriteria.locator('summary').click();
  await expect(page.locator('#criterionCaiMode')).toBeVisible();
  await page.locator('#sequenceInput').fill('>one\nATGTCCAAG\n>two\nATGTCCAAG');

  await expect(page.locator('#validationWarning')).toContainText('Multiple FASTA records detected');
  await expect(page.locator('#optimizeBtn')).toBeDisabled();
});

test('offers host selection while marking BY-2 experimental', async ({ page }) => {
  await openApp(page);

  await expect(page.locator('#hostSelect')).toHaveValue('nbenthamiana');
  await expect(page.locator('#hostSelect option[value="by2"]')).toContainText('experimental');
  await expect(page).toHaveTitle('FactorForge | N. benthamiana CDS Design');
  await expect(page.locator('input[name="objective"][value="feasibility_best"]')).toBeEnabled();
  await page.locator('#hostSelect').selectOption('by2');
  await expect(page.locator('#appliedPolicySummary')).toContainText('Tobacco BY-2');
});

test('renders experimental dual comparison and paged codon alignment', async ({ page }) => {
  await page.route('**/api/optimize', async route => {
    if (route.request().method() === 'GET') {
      await route.fulfill({ status: 200, contentType: 'application/json', body: JSON.stringify({
        capabilities: {
          execution_modes: ['profile', 'slm', 'dual_compare'],
          ml_preview: { available: true, status: 'experimental', trained_model_loaded: false },
          db_save: { available: false }
        },
        validation_checks: []
      }) });
      return;
    }
    const alignment = Array.from({ length: SAMPLE_PROTEIN.length }, (_, index) => ({
      position: index + 1,
      amino_acid: SAMPLE_PROTEIN[index],
      rule_codon: index === 3 ? 'CTT' : 'GAA',
      ml_codon: index === 3 ? 'TTA' : 'GAA',
      is_different: index === 3,
      rule_frequency: 0.25,
      ml_frequency: 0.15,
      rule_gc_bases: 1,
      ml_gc_bases: 0
    }));
    await route.fulfill({ status: 200, contentType: 'application/json', body: JSON.stringify({
      success: true,
      mode: 'dual_compare',
      experimental: true,
      preview_notice: 'ML update in progress: constrained-decoding research scaffold.',
      optimized_sequence: MOCK_DNA,
      original_length: SAMPLE_PROTEIN.length,
      optimized_length: MOCK_DNA.length,
      metrics: { cai: 0.91, gc_percent: 45, polya_signals: 0, length: MOCK_DNA.length },
      profile: 'balanced',
      host_profile: 'nbenthamiana',
      validation: { input_type: 'protein', polya: 'PASS', moclo: 'PASS', gc: 'PASS' },
      comparison: {
        rule: { cds: MOCK_DNA, metrics: { cai: 0.91, gc_percent: 45, type_iis_clean: true } },
        ml: { cds: MOCK_DNA, metrics: { cai: 0.88, gc_percent: 43, type_iis_clean: true } },
        codon_concordance_percent: 94.7,
        nt_identity_percent: 98.2,
        alignment,
        provenance: { canonical_sha256: 'a'.repeat(64), db_save_status: 'not_checked', audit_status: 'not_checked', leakage_check_status: 'not_checked' }
      }
    }) });
  });
  await openApp(page);
  await page.locator('input[name="engineMode"][value="dual_compare"]').check();
  await page.locator('#sequenceInput').fill(SAMPLE_PROTEIN);
  await page.locator('#optimizeBtn').click();

  await expect(page.locator('#comparisonDashboard')).toBeVisible();
  await expect(page.locator('#comparisonNotice')).toContainText('ML update in progress');
  await expect(page.locator('#comparisonMatrixBody')).toContainText('94.7% match');
  await expect(page.locator('#codonAlignmentViewer')).toContainText('▲');
  await expect(page.locator('#provenanceBadges')).toContainText('not_checked');
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
        profile: 'feasibility_best',
        host_profile: 'nbenthamiana',
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
  await page.locator('#optimizeBtn').click();

  await expect(page.locator('#resultsPanel')).toBeVisible();
  await expect(page.locator('#resultsPanel').getByRole('heading', { name: /Design Review/ })).toBeVisible();
  await expect.poll(() => requestBody).toMatchObject({
    sequence: SAMPLE_PROTEIN,
    host: 'nbenthamiana',
    host_profile: 'nbenthamiana',
    objective: 'feasibility_best',
    acceptance_criteria: expect.objectContaining({
      cai: expect.objectContaining({ mode: 'preferred' })
    })
  });
  await expect(page.locator('#hostProfileValue')).toContainText('nbenthamiana');
  await expect(page.locator('#optimizedSequence')).toContainText(MOCK_DNA.slice(0, 20));
});

test('optional seed and Type IIS presets are merged into the optimization payload', async ({ page }) => {
  let requestBody;
  await page.route('**/api/optimize', async route => {
    requestBody = route.request().postDataJSON();
    await route.fulfill({
      status: 200,
      contentType: 'application/json',
      body: JSON.stringify(reviewResponse({
        seed: 42,
        custom_restriction_sites: {
          requested: [{ name: 'SapI', sequence: 'GAAGAGC' }],
          detected: [],
          removed: [],
          unresolved: [],
        },
      }))
    });
  });
  await openApp(page);

  await page.locator('#sequenceInput').fill(SAMPLE_PROTEIN);
  const advancedSettings = page.locator('#advancedSettings');
  await expect(advancedSettings).not.toHaveAttribute('open', '');
  await advancedSettings.locator('summary').click();
  await page.locator('#optimizationSeed').fill('42');
  await page.locator('#customRestrictionSites').fill('SapI:GAAGAGC');
  await page.locator('input[name="typeIisEnzyme"][value="SapI"]').check();
  await page.locator('#optimizeBtn').click();

  await expect.poll(() => requestBody).toMatchObject({ seed: 42 });
  expect(requestBody.custom_restriction_sites.filter(site => site.name === 'SapI')).toHaveLength(1);
  await expect(page.locator('#customRestrictionResults')).toHaveClass(/domestication-attempted/);
  await expect(page.locator('#mfeWarningBanner')).toContainText('missing_dependency');
  await expect(page.locator('#gcTargetRange')).toHaveText('Target: 40.0–47.0%');
  await expect(page.locator('#resultsReportBody')).toContainText('seed=42');
  await expect(page.locator('#resultsReportBody')).toContainText('Type IIS');
  await expect(page.locator('#resultsReportBody')).toContainText('PASS');
  await expect(page.locator('#resultsReportBody')).toContainText('SapI');

  await page.locator('#resultsReport > summary').click();
  const [download] = await Promise.all([
    page.waitForEvent('download'),
    page.locator('#downloadResultsReportBtn').click(),
  ]);
  expect(download.suggestedFilename()).toBe('factorforge_design_review_ff-result-272.html');
  const downloadStream = await download.createReadStream();
  const chunks = [];
  for await (const chunk of downloadStream) chunks.push(chunk);
  const downloadedHtml = Buffer.concat(chunks).toString('utf-8');
  expect(downloadedHtml).toContain('<dt>Seed</dt><dd>42</dd>');
  expect(downloadedHtml).toContain('Researcher Decision Report');
  expect(downloadedHtml).toContain('Review priorities and next actions');
  expect(downloadedHtml).toContain('sha256:params-272');
});

test('results distinguish no domestication and hide the MFE warning when computed', async ({ page }) => {
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
        seed: null,
        metrics: {
          cai: 0.91,
          gc_percent: 45.0,
          polya_signals: 0,
          length: MOCK_DNA.length,
          mfe_status: 'computed',
          mfe_status_reason: null,
          requested_gc_min_percent: 40,
          requested_gc_max_percent: 47
        },
        profile: 'feasibility_best',
        host_profile: 'nbenthamiana',
        validation: { input_type: 'protein', polya: 'PASS', moclo: 'PASS', gc: 'PASS' }
      })
    });
  });
  await openApp(page);

  await page.locator('#sequenceInput').fill(SAMPLE_PROTEIN);
  await page.locator('#optimizeBtn').click();

  await expect.poll(() => requestBody && Object.hasOwn(requestBody, 'seed')).toBe(false);
  await expect(page.locator('#customRestrictionResults')).toHaveClass(/domestication-not-attempted/);
  await expect(page.locator('#customRestrictionResults')).toContainText('Domestication not attempted');
  await expect(page.locator('#mfeWarningBanner')).toBeHidden();
  await expect(page.locator('#resultsReportBody')).toContainText('seed not specified');
  await expect(page.locator('#resultsReportBody')).toContainText('MFE');
  await expect(page.locator('#resultsReportBody')).toContainText('Computed');
});

test('report treats the API decision and policy snapshot as authoritative', async ({ page }) => {
  const response = reviewResponse();
  response.metrics.cai = 0.75;
  response.qc_decision_matrix = [
    { criterion: 'cai', mode: 'preferred', observed: 0.75, threshold: 0.7, result: 'PASS' },
  ];
  response.acceptance_criteria_snapshot = { cai: { mode: 'preferred', minimum: 0.7 } };
  await mockOptimization(page, response);
  await openApp(page);

  await page.locator('#sequenceInput').fill(SAMPLE_PROTEIN);
  await page.locator('#optimizeBtn').click();

  const report = page.locator('#resultsReportBody');
  await expect(report).toContainText('Researcher decision brief');
  await expect(report).toContainText('PASS');
  await expect(report).toContainText('0.75');
  await expect(report).toContainText('0.7');
  await expect(report).not.toContainText('0.800 minimum');
});

test('report distinguishes preferred warnings and unavailable computation', async ({ page }) => {
  const response = reviewResponse({
    automated_decision: 'CONDITIONAL_PASS',
    decision_summary: { required_failure_count: 0, preferred_warning_count: 1, explanation: 'overall_gc requires review.' },
    qc_decision_matrix: [
      { criterion: 'overall_gc', mode: 'preferred', observed: 35.9, threshold: '40-47%', result: 'WARN' },
    ],
  });
  response.metrics.gc_percent = 35.9;
  await mockOptimization(page, response);
  await openApp(page);

  await page.locator('#sequenceInput').fill(SAMPLE_PROTEIN);
  await page.locator('#optimizeBtn').click();

  const report = page.locator('#resultsReportBody');
  await expect(report).toContainText('CONDITIONAL PASS');
  await expect(report.locator('[data-report-status="WARNING"]').first()).toContainText('WARNING');
  await expect(report.locator('[data-report-status="NOT_COMPUTED"]').first()).toContainText('NOT COMPUTED');
  await expect(report).toContainText('missing_dependency');
});

test('report preserves a required failure instead of softening it', async ({ page }) => {
  const response = reviewResponse({
    automated_decision: 'FAIL',
    decision_summary: { required_failure_count: 1, preferred_warning_count: 0, explanation: 'type_iis requires review.' },
    qc_decision_matrix: [
      { criterion: 'type_iis', mode: 'required', observed: 1, threshold: 0, result: 'FAIL' },
    ],
    acceptance_evaluation: {
      optimized: { criteria: [], details: { type_iis_sites: [{ enzyme: 'BsaI', site: 'GGTCTC', start: 12 }] } },
    },
  });
  await mockOptimization(page, response);
  await openApp(page);

  await page.locator('#sequenceInput').fill(SAMPLE_PROTEIN);
  await page.locator('#optimizeBtn').click();

  const report = page.locator('#resultsReportBody');
  await expect(report).toContainText('1 required fail');
  await expect(report.locator('[data-report-status="FAIL"]').first()).toContainText('FAIL');
  await expect(report).toContainText('BsaI at nt 13');
  await expect(report).toContainText('redesign or explicitly resolve every required site');
});

test('report exposes requested and applied Type IIS settings when they differ', async ({ page }) => {
  const response = reviewResponse({
    custom_restriction_sites: {
      requested: [{ name: 'BsaI', sequence: 'GGTCTC' }],
      detected: [],
      removed: [],
      unresolved: [],
    },
    acceptance_criteria_snapshot: {
      type_iis: { mode: 'required', enzymes: ['BsaI', 'BsmBI/Esp3I', 'SapI'] },
    },
  });
  await mockOptimization(page, response);
  await openApp(page);

  await page.locator('#sequenceInput').fill(SAMPLE_PROTEIN);
  await page.locator('#optimizeBtn').click();

  const report = page.locator('#resultsReportBody');
  await expect(report).toContainText('Type IIS settings mismatch');
  await expect(report).toContainText('Requested BsaI; applied BsaI, BsmBI/Esp3I, SapI.');
  await expect(report).toContainText('Mismatch');
});

test('evidence JSON matches the report and excludes raw sequences', async ({ page }) => {
  await mockOptimization(page, reviewResponse());
  await openApp(page);
  await page.locator('#sequenceInput').fill(SAMPLE_PROTEIN);
  await page.locator('#optimizeBtn').click();
  await page.locator('#resultsReport > summary').click();

  const [download] = await Promise.all([
    page.waitForEvent('download'),
    page.locator('#downloadEvidenceRecordBtn').click(),
  ]);
  expect(download.suggestedFilename()).toBe('factorforge_design_evidence_ff-result-272.json');
  const stream = await download.createReadStream();
  const chunks = [];
  for await (const chunk of stream) chunks.push(chunk);
  const text = Buffer.concat(chunks).toString('utf-8');
  const evidence = JSON.parse(text);
  expect(evidence.report_schema_version).toBe('1.0');
  expect(evidence.disposition.automated_decision).toBe('PASS');
  expect(evidence.provenance.parameter_hash).toBe('sha256:params-272');
  expect(evidence.artifacts).toBeUndefined();
  expect(text).not.toContain(MOCK_DNA);
  expect(text).not.toContain(SAMPLE_PROTEIN);
});

test('current history preserves report provenance without storing the raw input', async ({ page }) => {
  const response = reviewResponse({ input_type: 'cds', original_length: MOCK_DNA.length, validation: { input_type: 'cds', polya: 'PASS', moclo: 'PASS', gc: 'PASS' } });
  await mockOptimization(page, response);
  await openApp(page);
  await page.locator('#sequenceInput').fill(MOCK_DNA);
  await page.locator('#optimizeBtn').click();

  const stored = await page.evaluate(() => JSON.parse(localStorage.getItem('factorforge_history')));
  expect(stored.schemaVersion).toBe(3);
  expect(stored.items[0].inputSequence).toBeUndefined();
  expect(stored.items[0].resultSnapshot.provenance.parameter_hash).toBe('sha256:params-272');

  await page.locator('#historyList > div').first().click();
  await expect(page.locator('#sequenceInput')).toHaveValue('');
  await expect(page.locator('#resultsReportBody')).toContainText('sha256:params-272');
  await expect(page.locator('#resultsReportBody')).toContainText('nucleotide changes Not recorded');
});

test('legacy history remains readable and marks missing report fields', async ({ page }) => {
  await page.addInitScript(({ dna }) => {
    localStorage.setItem('factorforge_history', JSON.stringify({
      schemaVersion: 2,
      items: [{ id: 7, timestamp: 'legacy', inputLen: dna.length, profile: 'balanced', host: 'nbenthamiana', cai: 0.9, gc: 45, sequence: dna, inputSequence: dna }],
    }));
  }, { dna: MOCK_DNA });
  await openApp(page);

  await page.locator('#historyList > div').first().click();
  await expect(page.locator('#sequenceInput')).toHaveValue(MOCK_DNA);
  await expect(page.locator('#resultsReportBody')).toContainText('NOT AVAILABLE');
  await expect(page.locator('#resultsReportBody')).toContainText('Not recorded');
});

test('report remains usable in dark mode at a 390px viewport', async ({ page }) => {
  await page.setViewportSize({ width: 390, height: 844 });
  await mockOptimization(page, reviewResponse());
  await openApp(page);
  await page.locator('#themeToggle').click();
  await page.locator('#sequenceInput').fill(SAMPLE_PROTEIN);
  await page.locator('#optimizeBtn').click();
  await page.locator('#resultsReport > summary').click();

  await expect(page.locator('html')).toHaveClass(/dark/);
  await expect(page.locator('#design-review-report-title')).toBeVisible();
  await expect(page.locator('#downloadEvidenceRecordBtn')).toBeVisible();

  const downloadPromise = page.waitForEvent('download');
  await page.locator('#downloadResultsReportBtn').click();
  const download = await downloadPromise;
  const stream = await download.createReadStream();
  let reportHtml = '';
  for await (const chunk of stream) reportHtml += chunk.toString();
  expect(reportHtml).toContain('<html lang="en" data-theme="dark">');
  expect(reportHtml).toContain('[data-theme="dark"] body');
});

test('desktop design columns scroll with the page instead of sticking independently', async ({ page }) => {
  await page.setViewportSize({ width: 1440, height: 900 });
  await openApp(page);

  await expect(page.locator('#designBriefPanel')).toHaveCSS('position', 'static');
  await expect(page.locator('#resultsPanel')).toHaveCSS('position', 'static');
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
