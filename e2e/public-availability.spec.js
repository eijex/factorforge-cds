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

test('4DP displays sidecar positions before review and clears review after a replacement', async ({ page }) => {
  await page.route('**/api/optimize', route => route.fulfill({
    json: { capabilities: { local_4dp: { available: true, fixtures: false } } }
  }));
  await page.goto('/');
  await page.locator('input[value="paired_4dp"]').check();
  const review = page.locator('#review4dpPositions');
  const upload = page.locator('#pauseSidecarFile');
  const load = data => upload.setInputFiles({
    name: 'test-sidecar.json', mimeType: 'application/json',
    buffer: Buffer.from(JSON.stringify(data))
  });
  await expect(review).toBeDisabled();
  await load({
    primary_chain: 'A', tunnel_offset_applied: 0,
    chains: { A: { length: 50, pause_candidates: [
      { within_structure_chain: true, tunnel_pause_position: 38 },
      { within_structure_chain: false, tunnel_pause_position: 55 }
    ] } }
  });
  await expect(page.locator('#sidecarPositionReview')).toHaveText('Chain A · Offset 0 residues · Model positions: 38');
  await review.check();
  await load({});
  await expect(review).not.toBeChecked();
  await expect(review).toBeDisabled();
  await expect(page.locator('#loadSampleSidecarBtn')).toBeHidden();
  await expect(page.locator('#allowMockSidecarLabel')).toBeHidden();
});

test('3DP requires structure upload cleanly without dead sample references', async ({ page }) => {
  await page.route('**/api/optimize', route => route.fulfill({
    json: { capabilities: { direct_structure: { available: true, status: 'experimental' }, local_4dp: { available: true, fixtures: false } } }
  }));
  await page.goto('/');
  await page.locator('input[value="paired_3dp"]').check();
  await expect(page.locator('#structureWebInputs summary')).toContainText('Required: PDB / mmCIF structure upload and position review');
  await page.locator('#sequenceInput').fill('MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG');
  await page.locator('#optimizeBtn').click();
  await expect(page.locator('body')).toContainText('Upload an exact-match PDB or mmCIF structure file before generating 3DP.');
  await expect(page.locator('body')).not.toContainText('Adalimumab (Humira) 3DP Sample');
});

test('drag and drop loads FASTA, 3DP, and 4DP files cleanly without navigating away', async ({ page }) => {
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

  // Mismatched .pdb dropped onto Sequence Input shows guidance toast and keeps URL intact
  await dispatchDrop('#sequenceUploadZone', '1UBQ.pdb', 'chemical/x-pdb', 'ATOM      1  CA  MET A   1\n');
  await expect(page.locator('body')).toContainText('PDB/mmCIF structure files belong in 3DP');
  expect(page.url()).toContain('http://127.0.0.1:8010');

  // Valid .fasta dropped onto Sequence Input loads sequence
  await dispatchDrop('#sequenceUploadZone', 'test.fasta', 'text/plain', '>test\nMASTWQHLEKDPVNRGYF\n');
  await expect(page.locator('#sequenceInput')).toHaveValue(/MASTWQHLEKDPVNRGYF/);

  // 3DP .pdb dropped onto structure drop zone attaches to #structureWebFile
  await page.locator('input[value="paired_3dp"]').check();
  await dispatchDrop('#structureWebDropZone', '1UBQ.pdb', 'chemical/x-pdb', 'ATOM      1  CA  MET A   1       1.000   0.000   0.000  1.00 20.00           C  \n');
  const attached3dpName = await page.locator('#structureWebFile').evaluate(el => el.files?.[0]?.name || '');
  expect(attached3dpName).toBe('1UBQ.pdb');

  // 4DP .json dropped onto sidecar drop zone populates #sidecarPositionReview
  await page.locator('input[value="paired_4dp"]').check();
  await dispatchDrop('#sidecarDropZone', 'sidecar.json', 'application/json', JSON.stringify({
    primary_chain: 'A', tunnel_offset_applied: 35,
    chains: { A: { length: 76, pause_candidates: [{ within_structure_chain: true, tunnel_pause_position: 76 }] } }
  }));
  await expect(page.locator('#sidecarPositionReview')).toHaveText('Chain A · Offset 35 residues · Model positions: 76');
});
