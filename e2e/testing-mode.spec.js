const {test,expect}=require('@playwright/test');

// Retained custom SOP code is tested behind its disabled production feature flag.
test.beforeEach(async ({page}) => {
    const source = require('node:fs').readFileSync(require('node:path').join(__dirname, '../web/js/app.js'), 'utf8');
    await page.route('**/js/app.js', route => route.fulfill({contentType:'text/javascript', body:source.replace('const ENABLE_CUSTOM_SOP = false;', 'const ENABLE_CUSTOM_SOP = true;')}));
});
test('watermark remains disabled after importing an enabled profile', async({page})=>{
  let payload;
  await page.route('**/api/optimize',route=>{
    if(route.request().method()==='POST') {payload=route.request().postDataJSON();return route.fulfill({status:400,json:{error:'Synthetic capture'}});}
    return route.fulfill({json:{supported_objectives:['dp_v2_1_1'],capabilities:{},validation_checks:[]}});
  });
  await page.goto('/');
  await expect(page.locator('#researchTestingStatus')).toContainText('Both modes use this form');
  await page.locator('#experimentalSettings > summary').click();
  await expect(page.locator('#toggleWatermark')).toBeDisabled();
  await expect(page.locator('#toggleWatermark')).not.toBeChecked();
  // Use the currently rendered profile export contract for a valid upload.
  const custom={$schema:'factorforge-sop-v1',profile_id:'watermark_test',sop_name:'Synthetic test',version:'1.0.0',author:'Synthetic',status:'LOCAL',default_enforcement:'IGNORE',unknown_rule_policy:'ERROR',rules:{},workflow:{watermark:true,design_method:'feasibility_best',comparison_methods:[],sequence_requirements:{type_iis_enzymes:['BsaI'],custom_restriction_sites:[],reproducibility_seed:null,kozak_optimization:false,suppress_tpa:false,moclo_template:false},review_policy:{cai:'preferred',overall_gc:'required',local_gc:'preferred',type_iis:'required',repeats:'preferred',homopolymers:'preferred',forbidden_motifs:'preferred'}}};
  await page.locator('#sopFileUpload').setInputFiles({name:'synthetic.json',mimeType:'application/json',buffer:Buffer.from(JSON.stringify(custom))});
  await expect(page.locator('#sopProfileName')).toHaveText('Synthetic test');
  await expect(page.locator('#toggleWatermark')).toBeDisabled();
  await expect(page.locator('#toggleWatermark')).not.toBeChecked();
  await page.locator('#sequenceInput').fill('MASTWQHLEKDPVNRSGYF');
  await expect(page.locator('#optimizeBtn')).toBeEnabled();
  await page.locator('#optimizeBtn').click();
  await expect.poll(()=>payload).toBeTruthy();
  expect(payload.sop_profile.workflow.watermark).toBe(false);
  await page.reload();
  await expect(page.locator('#toggleWatermark')).toBeDisabled();
  await expect(page.locator('#toggleWatermark')).not.toBeChecked();
});

test('3DP auto-resolves or accepts optional structure input and submits to the same API without opening ResearchOS', async ({page}) => {
  let submitted=null;
  await page.route('**/api/optimize',route=>{
    if(route.request().method()==='GET') return route.fulfill({json:{capabilities:{direct_structure:{available:true}},validation_checks:[]}});
    const body=route.request().postDataJSON();
    if(body?.action==='auto_resolve_structure') return route.fulfill({status:200,json:{success:true,action:'auto_resolve_structure',resolved:{content:'ATOM      1  CA  MET A   1       0.000   0.000   0.000  1.00 85.00           C  \nEND\n',format:'pdb',chain:'A',model:1,positions:[3,4,5],positions_reviewed:true,source_label:'Synthetic Auto-Resolved',ca_trace:[]}}});
    submitted=body;
    return route.fulfill({status:400,json:{error:'Synthetic structure rejection'}});
  });
  await page.goto('/');
  await page.evaluate(()=>{window.open=()=>{throw new Error('Unexpected ResearchOS navigation');};});
  await page.locator('input[name="researchDesignMode"][value="paired_3dp"]').check();
  await expect(page.locator('#structureWebInputs')).toBeVisible();
  await expect(page.locator('#appliedPolicySummary')).toContainText('3DP selected');
  await expect(page.locator('#optimizeBtn')).toBeDisabled();
  await page.locator('#sequenceInput').fill('MASTWQHLEKDPVNRSGYF');
  await page.locator('#structureWebFile').setInputFiles({name:'synthetic.pdb',mimeType:'text/plain',buffer:Buffer.from('synthetic test content')});
  await page.locator('#structureWebPositions').fill('3,4,5');
  await page.locator('#structureWebReviewed').check();
  await page.locator('#optimizeBtn').click();
  await expect.poll(()=>submitted).toBeTruthy();
  expect(submitted.design_mode).toBe('paired_3dp');
  expect(submitted.structure_context.positions).toEqual([3,4,5]);
  expect(submitted.kozak).toBe(false);
  expect(submitted.sop_profile.workflow.watermark).toBe(false);
  await expect(page.locator('#emptyState')).not.toHaveClass(/\bhidden\b/);
  await expect(page.locator('#candidateComparisonBody')).toBeEmpty();
  await page.locator('input[name="researchDesignMode"][value="baseline_1dp"]').check();
  await expect(page.locator('#structureWebInputs')).toBeHidden();
  await expect(page.locator('#appliedPolicySummary')).toContainText('1DP selected');
});

test('3DP settings export and import restore controls without restoring structure or approval', async({page})=>{
  const fs=require('node:fs');
  await page.route('**/api/optimize',r=>r.fulfill({json:{capabilities:{direct_structure:{available:true}},validation_checks:[]}}));
  await page.goto('/');
  await page.locator('input[name="researchDesignMode"][value="paired_3dp"]').check();
  await page.locator('#structureWebChain').fill('B');
  await page.locator('#structureWebPositions').fill('5,6');
  await page.locator('#structureWebPenalty').fill('3');
  await page.locator('#structureWebReviewed').check();
  await page.locator('#profileManagement > summary').click();
  const downloaded=page.waitForEvent('download');
  await page.locator('#downloadSopTemplate').click();
  const file=await downloaded;const text=fs.readFileSync(await file.path(),'utf8');
  expect(text).toContain('design_mode: paired_3dp');
  expect(text).toContain('positions: [5, 6]');
  expect(text).not.toContain('positions_reviewed');
  await page.locator('input[name="researchDesignMode"][value="baseline_1dp"]').check();
  await page.locator('#sopFileUpload').setInputFiles({name:'roundtrip.yaml',mimeType:'application/yaml',buffer:Buffer.from(text)});
  await expect(page.locator('input[name="researchDesignMode"][value="paired_3dp"]')).toBeChecked();
  await expect(page.locator('#structureWebInputs')).toBeVisible();
  await expect(page.locator('#structureWebChain')).toHaveValue('B');
  await expect(page.locator('#structureWebPositions')).toHaveValue('5,6');
  await expect(page.locator('#structureWebPenalty')).toHaveValue('3');
  await expect(page.locator('#structureWebReviewed')).not.toBeChecked();
  await expect(page.locator('#structureWebFile')).toHaveValue('');
});
