const { test, expect } = require('@playwright/test');

test('public 1DP/3DP/4DP availability and custom SOP collaboration terms remain consistent', async ({ page }, testInfo) => {
  await page.route('**/api/optimize', route => route.fulfill({
    json: { capabilities: { direct_structure: { available: true, status: 'experimental' }, local_4dp: { available: true, fixtures: false } } }
  }));
  await page.goto('/');
  await expect(page.locator('input[value="baseline_1dp"]')).toBeChecked();
  await expect(page.locator('input[value="paired_3dp"]')).toBeEnabled();
  await expect(page.locator('input[value="paired_4dp"]')).toBeEnabled();
  await expect(page.locator('#policyScopeNotice')).toContainText('not an approved laboratory SOP');
  await expect(page.locator('#policyScopeNotice')).toContainText('Agree on sharing');
  await expect(page.locator('body')).toContainText('eijex.lab@gmail.com');
  await expect(page.locator('body')).toContainText('GNU AGPL-3.0 applies to the FactorForge software');
  await page.screenshot({ path: testInfo.outputPath('public-policy.png'), fullPage: true });
});

test('explicit API pause disables 3DP and 4DP cleanly', async ({ page }) => {
  await page.route('**/api/optimize', route => route.fulfill({
    json: { capabilities: { direct_structure: { available: false, status: 'paused' }, local_4dp: { available: false, fixtures: false } } }
  }));
  await page.goto('/');
  await expect(page.locator('input[value="baseline_1dp"]')).toBeChecked();
  await expect(page.locator('input[value="paired_3dp"]')).toBeDisabled();
  await expect(page.locator('input[value="paired_4dp"]')).toBeDisabled();
  await page.evaluate(() => {
    const input = document.querySelector('input[value="paired_4dp"]');
    input.disabled = false;
    input.checked = true;
    input.dispatchEvent(new Event('change'));
  });
  await expect(page.locator('input[value="baseline_1dp"]')).toBeChecked();
  await expect(page.locator('#structure4dpInputs')).toBeHidden();
});

test('4DP displays verified sidecar positions before review and clears review after a replacement', async ({ page }) => {
  await page.route('**/api/optimize', route => route.fulfill({
    json: { capabilities: { local_4dp: { available: true, fixtures: false } } }
  }));
  await page.goto('/');
  const seq = 'MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLG';
  await page.locator('#sequenceInput').fill(seq);
  await page.locator('input[value="paired_4dp"]').check();
  const review = page.locator('#review4dpPositions');
  const ack = page.locator('[data-role="ack4dpHypothesis"]');
  const upload = page.locator('#pauseSidecarFile');
  const load = data => upload.setInputFiles({
    name: 'test-sidecar.json', mimeType: 'application/json',
    buffer: Buffer.from(JSON.stringify(data))
  });
  await expect(review).toBeDisabled();
  await expect(ack).toBeDisabled();

  // Mismatched/incomplete sidecar remains blocked
  await load({
    primary_chain: 'A', tunnel_offset_applied: 0,
    chains: { A: { length: 50, pause_candidates: [{ within_structure_chain: true, tunnel_pause_position: 38 }] } }
  });
  await expect(review).toBeDisabled();
  await expect(page.locator('#sidecarPositionReview')).toContainText('Blocked:');

  // Complete integrity-passing sidecar unlocks review & hypothesis acknowledgment
  await load({
    status: 'AVAILABLE',
    coordinate_system: 'mmCIF_label_seq_id_1based',
    source_type: 'ESM_ATLAS_REMOTE',
    primary_chain: 'A',
    tunnel_offset_applied: 0,
    chains: {
      A: {
        sequence: seq,
        length: seq.length,
        residues: seq.split('').map((aa, idx) => ({ seq_id: idx + 1, aa, plddt: 88.5 })),
        pause_candidates: [
          { within_structure_chain: true, tunnel_pause_position: 38 },
          { within_structure_chain: false, tunnel_pause_position: 55 }
        ]
      }
    }
  });
  await expect(page.locator('#sidecarPositionReview')).toHaveText('Chain A · Offset 0 residues · Verified positions: 38');
  await expect(review).toBeEnabled();
  await expect(ack).toBeEnabled();
  await review.check();
  await ack.check();
  await load({});
  await expect(review).not.toBeChecked();
  await expect(review).toBeDisabled();
  await expect(ack).not.toBeChecked();
  await expect(ack).toBeDisabled();
  await expect(page.locator('#loadSampleSidecarBtn')).toBeHidden();
  await expect(page.locator('#allowMockSidecarLabel')).toBeHidden();
});

test('3DP requires structure upload cleanly without dead sample references', async ({ page }) => {
  await page.route('**/api/optimize', route => route.fulfill({
    json: { capabilities: { direct_structure: { available: true, status: 'experimental' }, local_4dp: { available: true, fixtures: false } } }
  }));
  await page.goto('/');
  await page.locator('input[value="paired_3dp"]').check();
  await expect(page.locator('#structureWebInputs summary')).toContainText('Required: PDB / mmCIF structure upload and chain/position review');
  await page.locator('#sequenceInput').fill('MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG');
  await page.locator('#optimizeBtn').click();
  await expect(page.locator('body')).toContainText('Upload an exact-match PDB or mmCIF structure file before generating 3DP.');
  await expect(page.locator('body')).not.toContainText('Adalimumab (Humira) 3DP Sample');
});

test('unified drag and drop loads FASTA, 3DP, and 4DP files cleanly without navigating away', async ({ page }) => {
  await page.route('**/api/optimize', route => route.fulfill({
    json: { capabilities: { direct_structure: { available: true, status: 'experimental' }, local_4dp: { available: true, fixtures: false } } }
  }));
  await page.goto('/');

  const dispatchDrop = async (selector, name, mimeType, content) => {
    await page.evaluate(async ({ sel, fileName, type, text }) => {
      const target = document.querySelector(sel);
      const dt = new DataTransfer();
      dt.items.add(new File([text], fileName, { type }));
      target.dispatchEvent(new DragEvent('dragover', { bubbles: true, cancelable: true, dataTransfer: dt }));
      target.dispatchEvent(new DragEvent('drop', { bubbles: true, cancelable: true, dataTransfer: dt }));
    }, { sel: selector, fileName: name, type: mimeType, text: content });
  };

  // Valid .fasta dropped onto unified Sequence Input loads sequence
  await dispatchDrop('#sequenceUploadZone', 'test.fasta', 'text/plain', '>test\nMASTWQHLEKDPVNRGYF\n');
  await expect(page.locator('#sequenceInput')).toHaveValue(/MASTWQHLEKDPVNRGYF/);

  // .pdb dropped onto unified Sequence Input routes directly to 3DP without navigating away
  await dispatchDrop('#sequenceUploadZone', '1UBQ.pdb', 'chemical/x-pdb', 'ATOM      1  CA  MET A   1       1.000   0.000   0.000  1.00 20.00           C  \n');
  const attached3dpName = await page.locator('#structureWebFile').evaluate(el => el.files?.[0]?.name || '');
  expect(attached3dpName).toBe('1UBQ.pdb');
  expect(page.url()).toContain('http://127.0.0.1:8010');

  // 4DP .json dropped onto unified zone validates sidecar against active sequence
  const seq = 'MASTWQHLEKDPVNRGYF';
  await dispatchDrop('#sequenceUploadZone', 'sidecar.json', 'application/json', JSON.stringify({
    status: 'AVAILABLE',
    coordinate_system: 'mmCIF_label_seq_id_1based',
    source_type: 'ESM_ATLAS_REMOTE',
    primary_chain: 'A',
    tunnel_offset_applied: 10,
    chains: {
      A: {
        sequence: seq,
        length: seq.length,
        residues: seq.split('').map((aa, idx) => ({ seq_id: idx + 1, aa, plddt: 90.0 })),
        pause_candidates: [{ within_structure_chain: true, tunnel_pause_position: 12 }]
      }
    }
  }));
  await expect(page.locator('#sidecarPositionReview')).toHaveText('Chain A · Offset 10 residues · Verified positions: 12');
});

test('multi-FASTA, multi-chain PDB inspection, bounded Compare queue, and matrix labels work with synthetic fixtures', async ({ page }) => {
  const seqA = 'MASTWQHLEKDPVNRGYFMA';
  const seqB = 'MKTVRQERLKSIVRILERSK';
  let inFlight = 0;
  let maxInFlight = 0;

  await page.route('**/api/optimize', async route => {
    const req = route.request();
    if (req.method() === 'GET') {
      return route.fulfill({
        json: { capabilities: { direct_structure: { available: true, status: 'experimental' }, local_4dp: { available: true, fixtures: false } } }
      });
    }
    inFlight++;
    if (inFlight > maxInFlight) maxInFlight = inFlight;
    const body = req.postDataJSON();
    await new Promise(r => setTimeout(r, 40));
    inFlight--;

    if (body.action === 'auto_resolve_structure') {
      return route.fulfill({
        json: { resolved: { positions: [1, 2], source_label: 'Synthetic test fixture' } }
      });
    }

    const aaSeq = body.sequence;
    const mode = body.design_mode || '1dp';
    const codon12 = mode === 'paired_3dp' ? 'ATGGCTAGCACC' : 'ATGGCATCCACC';
    const fullDna = (codon12 + 'TGGCAACATCTTGAGAAGGACCCTGTGAACCGTGGTTACTTCATGGCT'.repeat(2)).slice(0, aaSeq.length * 3);
    return route.fulfill({
      json: {
        success: true,
        optimized_sequence: fullDna,
        original_length: aaSeq.length,
        optimized_length: fullDna.length,
        input_type: 'protein',
        metrics: {
          cai: mode === 'paired_3dp' ? 0.88 : 0.91,
          gc_percent: 44.5,
          polya_signals: 0,
          length: fullDna.length,
          mfe_kcal_mol: null,
          mfe_status: 'not_computed',
          mfe_used: false,
          requested_gc_min_percent: 40,
          requested_gc_max_percent: 47,
        },
        profile: mode === 'paired_3dp' ? 'balanced' : 'feasibility_best',
        host_profile: 'nbenthamiana',
        automated_decision: 'PASS',
        decision_summary: { required_failure_count: 0, preferred_warning_count: 0, explanation: 'All active acceptance criteria passed.' },
        validation: { input_type: 'protein', polya: 'PASS', moclo: 'PASS', gc: 'PASS' },
        constraint_report: { aa_identity: 1 },
        construct_id: `SYN-${mode.toUpperCase()}`,
        product_version: '3.6.4',
        engine_versions: { product: '3.6.4' },
        design_mode: mode,
        evidence_classification: 'UNVERIFIED_PUBLIC_CORPUS',
        biological_validation: 'NOT_AVAILABLE'
      }
    });
  });

  await page.goto('/');

  // Upload 2-record synthetic Multi-FASTA
  await page.locator('#sequenceInput').fill(`>ChainA_Syn\n${seqA}\n>ChainB_Syn\n${seqB}`);
  await page.locator('#sequenceInput').dispatchEvent('input');
  await expect(page.locator('#projectRecordsPanel')).toBeVisible();
  await expect(page.locator('#projectRecordsPanel')).toContainText('2 sequences detected');

  // Attach synthetic 2-chain PDB matching ChainA_Syn on Chain A and partial mismatch on Chain B
  const pdbLines = [
    ...seqA.split('').map((_, idx) => {
      const resn = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'];
      const map1to3 = { A:'ALA',R:'ARG',N:'ASN',D:'ASP',C:'CYS',Q:'GLN',E:'GLU',G:'GLY',H:'HIS',I:'ILE',L:'LEU',K:'LYS',M:'MET',F:'PHE',P:'PRO',S:'SER',T:'THR',W:'TRP',Y:'TYR',V:'VAL' };
      const aa3 = map1to3[_] || resn[0];
      return `ATOM  ${String(idx + 1).padStart(5)}  CA  ${aa3} A${String(idx + 1).padStart(4)}       1.000   0.000   0.000  1.00 85.00           C`;
    }),
    'TER',
    'ATOM     21  CA  MET B   1       2.000   0.000   0.000  1.00 70.00           C',
    'END'
  ].join('\n');
  await page.locator('#structureWebFile').setInputFiles({
    name: 'synthetic_h2.pdb', mimeType: 'chemical/x-pdb', buffer: Buffer.from(pdbLines)
  });
  await expect(page.locator('#structureChainInspector')).toBeVisible();
  await expect(page.locator('#structureChainInspector')).toContainText('Chain A');
  await expect(page.locator('#structureChainInspector')).toContainText('Chain B');

  // Switch to Compare mode, run, check bounded queue & truthful results matrix
  await page.locator('[data-strategy-mode="compare"]').click();
  await expect(page.locator('#compareMethodsBar')).toBeVisible();
  await page.locator('#structureWebPositions').fill('1,2');
  await page.locator('#structureWebReviewed').check();
  await page.locator('#optimizeBtn').click();

  await expect(page.locator('#resultsPanel')).toBeVisible();
  expect(maxInFlight).toBeLessThanOrEqual(2);
  await expect(page.locator('#projectChainTabs')).toContainText('ChainA_Syn');
  await expect(page.locator('#projectChainTabs')).toContainText('ChainB_Syn');
  await expect(page.locator('#matrixTableBody')).toContainText(/Cryptic splice prediction\s*NOT RUN/);
  await expect(page.locator('#matrixTableBody')).toContainText(/Biological validation\s*NOT AVAILABLE/);
  await expect(page.locator('#multiCandidateCodonDiffContainer')).toContainText(/2\s*·\s*A/);
  await expect(page.locator('#multiCandidateCodonDiffContainer')).toContainText('Synonymous codon difference');
  await expect(page.locator('#multiCandidateCodonDiffContainer')).toContainText('Recorded structure/4DP pause target');
});
