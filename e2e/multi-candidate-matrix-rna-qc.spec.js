const { test, expect } = require('@playwright/test');

test.beforeEach(async ({ page }) => {
    await page.route('**/api/optimize', route => route.fulfill({
        json: {
            capabilities: {
                direct_structure: { available: false },
                structure_recommendation: { available: false },
            },
        },
    }));
});

test('1DP / 3DP / 4DP matrix aligns RNA motif QC with each primary candidate (not candidates[0])', async ({ page }) => {
    await page.goto('/');
    await page.waitForFunction(() => typeof renderMultiCandidateMatrix === 'function');

    await page.evaluate(() => {
        document.getElementById('resultsContainer')?.classList.remove('hidden');
        state.projectResults = [
            {
                record: { id: 'rec_1', name: 'TestProt', effectiveType: 'protein', sequence: 'MASTNK' },
                resultsByMode: {
                    baseline_1dp: {
                        design_mode: 'baseline_1dp',
                        candidates: [
                            {
                                id: 'baseline_1dp',
                                dna_sequence: 'ATGGCTTCAACTAATAAATGA',
                                cai: 0.91,
                                gc_percent: 42.0,
                                scientific_review: {
                                    rna_risks: {
                                        summary_status: 'ADVISORY_FINDINGS',
                                        donor_count: 0,
                                        polya_count: 1,
                                    },
                                },
                            },
                        ],
                        recommended_candidate: {
                            id: 'baseline_1dp',
                            dna_sequence: 'ATGGCTTCAACTAATAAATGA',
                            cai: 0.91,
                            gc_percent: 42.0,
                            scientific_review: {
                                rna_risks: {
                                    summary_status: 'ADVISORY_FINDINGS',
                                    donor_count: 0,
                                    polya_count: 1,
                                },
                            },
                        },
                        scientific_review: {
                            rna_risks: { summary_status: 'ADVISORY_FINDINGS', donor_count: 0, polya_count: 1 },
                        },
                    },
                    paired_3dp: {
                        design_mode: 'paired_3dp',
                        // Intentionally place baseline_1dp in candidates[0] with distinct counts (2 splice, 1 polyA)
                        // while 3DP primary is candidates[1] with CLEAN_SCAN (0 splice, 0 polyA).
                        candidates: [
                            {
                                id: 'baseline_1dp',
                                dna_sequence: 'ATGGCTTCAACTAATAAATGA',
                                cai: 0.91,
                                gc_percent: 42.0,
                                scientific_review: {
                                    rna_risks: { summary_status: 'ADVISORY_FINDINGS', donor_count: 2, polya_count: 1 },
                                },
                            },
                            {
                                id: 'paired_3dp',
                                dna_sequence: 'ATGGCATCTACCAACAAGTGA',
                                cai: 0.89,
                                gc_percent: 45.0,
                                scientific_review: {
                                    rna_risks: { summary_status: 'CLEAN_SCAN', donor_count: 0, polya_count: 0 },
                                },
                            },
                        ],
                        recommended_candidate: {
                            id: 'paired_3dp',
                            dna_sequence: 'ATGGCATCTACCAACAAGTGA',
                            cai: 0.89,
                            gc_percent: 45.0,
                            scientific_review: {
                                rna_risks: { summary_status: 'CLEAN_SCAN', donor_count: 0, polya_count: 0 },
                            },
                        },
                        scientific_review: {
                            rna_risks: { summary_status: 'CLEAN_SCAN', donor_count: 0, polya_count: 0 },
                        },
                    },
                    paired_4dp: {
                        design_mode: 'paired_4dp',
                        // Test matchedCand lookup when the recommended candidate summary disagrees with its detailed scan
                        // and candidates[0] has 0/0 while candidates[1] (4DP) has 1 splice / 2 polyA.
                        candidates: [
                            {
                                id: 'baseline_1dp',
                                dna_sequence: 'ATGGCTTCAACTAATAAATGA',
                                cai: 0.91,
                                gc_percent: 42.0,
                                scientific_review: {
                                    rna_risks: { summary_status: 'CLEAN_SCAN', donor_count: 0, polya_count: 0 },
                                },
                            },
                            {
                                id: 'structure_4dp',
                                dna_sequence: 'ATGGCCAGCACAAACAAGTGA',
                                cai: 0.87,
                                gc_percent: 46.0,
                                scientific_review: {
                                    rna_risks: { summary_status: 'ADVISORY_FINDINGS', donor_count: 1, polya_count: 2 },
                                },
                            },
                        ],
                        recommended_candidate: {
                            id: 'structure_4dp',
                            dna_sequence: 'ATGGCCAGCACAAACAAGTGA',
                            scientific_review: { rna_risks: { summary_status: "CLEAN_SCAN", donor_count: 0, polya_count: 0 } },
                            cai: 0.87,
                            gc_percent: 46.0,
                        },
                    },
                },
            },
        ];
        state.activeResultRecordIndex = 0;
        renderMultiCandidateMatrix();
    });

    const row = label => page.locator('#matrixTableBody tr').filter({ hasText: label });

    // Columns: 0=Feature, 1=1DP, 2=3DP, 3=4DP
    const spliceCells = row('Splice-like motif scan').locator('td');
    await expect(spliceCells.nth(1)).toHaveText('RUN — 0 hits');
    await expect(spliceCells.nth(2)).toHaveText('RUN — 0 hits');
    await expect(spliceCells.nth(3)).toHaveText('RUN — 1 hit');
    await expect(row('Splice-like motif scan')).not.toContainText('PASS');

    const polyaCells = row('Poly-A motif scan').locator('td');
    await expect(polyaCells.nth(1)).toHaveText('RUN — 1 hit');
    await expect(polyaCells.nth(2)).toHaveText('RUN — 0 hits');
    await expect(polyaCells.nth(3)).toHaveText('RUN — 2 hits');
    await expect(row('Poly-A motif scan')).not.toContainText('PASS');

    // Separate prediction and wet-lab rows stay NOT RUN / NOT AVAILABLE
    for (const idx of [1, 2, 3]) {
        await expect(row('Cryptic splice prediction').locator('td').nth(idx)).toHaveText('NOT RUN');
        await expect(row('Biological validation').locator('td').nth(idx)).toHaveText('NOT AVAILABLE');
    }
});

test('strict RNA status allowlist rejects unknown statuses and avoids cross-candidate fallback', async ({ page }) => {
    await page.goto('/');
    await page.waitForFunction(() => typeof renderMultiCandidateMatrix === 'function');

    await page.evaluate(() => {
        document.getElementById('resultsContainer')?.classList.remove('hidden');
        state.projectResults = [];
        state.activeResultRecordIndex = 0;
        state.results = {
            input_type: 'protein',
            original_sequence: 'MASTNK',
            scientific_review: {
                rna_risks: { summary_status: 'CLEAN_SCAN', donor_count: 0, polya_count: 0 },
            },
            candidates: [
                {
                    id: 'cand_clean',
                    label: 'Candidate A',
                    dna_sequence: 'ATGGCATCTACCAACAAGTGA',
                    cai: 0.90,
                    gc_percent: 45.0,
                    scientific_review: {
                        rna_risks: { summary_status: 'CLEAN_SCAN', donor_count: 0, polya_count: 0 },
                    },
                },
                {
                    id: 'cand_unknown_status',
                    label: 'Candidate B',
                    dna_sequence: 'ATGGCATCCACTAACAAGTGA',
                    cai: 0.88,
                    gc_percent: 44.0,
                    scientific_review: {
                        // Unknown status must remain NOT RUN even if numeric counts are present
                        rna_risks: { summary_status: 'PENDING', donor_count: 0, polya_count: 0 },
                    },
                },
                {
                    id: 'cand_missing_review',
                    label: 'Candidate C',
                    dna_sequence: 'ATGGCCTCAACTAACAAGTGA',
                    cai: 0.86,
                    gc_percent: 43.0,
                },
            ],
        };
        renderMultiCandidateMatrix();
    });

    const row = label => page.locator('#matrixTableBody tr').filter({ hasText: label });
    const spliceCells = row('Splice-like motif scan').locator('td');
    const polyaCells = row('Poly-A motif scan').locator('td');

    await expect(spliceCells.nth(1)).toHaveText('RUN — 0 hits');
    await expect(polyaCells.nth(1)).toHaveText('RUN — 0 hits');

    await expect(spliceCells.nth(2)).toHaveText('NOT RUN');
    await expect(polyaCells.nth(2)).toHaveText('NOT RUN');

    await expect(spliceCells.nth(3)).toHaveText('NOT RUN');
    await expect(polyaCells.nth(3)).toHaveText('NOT RUN');
});
