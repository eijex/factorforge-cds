const {test,expect}=require('@playwright/test');
test.beforeEach(async({page})=>{
 await page.route('**/api/optimize',route=>route.fulfill({json:route.request().method()==='GET'?
 {capabilities:{direct_structure:{available:true},structure_recommendation:{available:true}}}:
 {success:true,recommendation:{status:'NOT_AVAILABLE',reason:'Unknown confidence source',recommended_positions:[],reviewed:false}}}));
 await page.goto('/');
 await page.waitForFunction(()=>apiCapabilities.structure_recommendation?.available);
});
test('unsupported structure observation never populates or approves positions',async({page})=>{
 await page.locator('input[name="researchDesignMode"][value="paired_3dp"]').check();
 await page.locator('#sequenceInput').fill('MASTWQHLEKDPVNRSGYF');
 await page.locator('#structureWebFile').setInputFiles({name:'synthetic.pdb',mimeType:'text/plain',buffer:Buffer.from('END')});
 await page.locator('#recommendStructurePositions').click();
 await expect(page.locator('#structureRecommendationStatus')).toContainText('NOT_AVAILABLE');
 await expect(page.locator('#acceptStructurePositions')).toBeDisabled();
 await expect(page.locator('#structureWebPositions')).toHaveValue('');
 await expect(page.locator('#structureWebReviewed')).not.toBeChecked();
});
test('missing observations display unavailable rather than zero or PASS',async({page})=>{
 const html=await page.evaluate(()=>scientificReviewHtml({scientific_review:{initiation_mfe:{status:'NOT_AVAILABLE',mfe_kcal_mol:null}}}));
 expect(html).toContain('Not measured');
 expect(html).toContain('donor N/A');
 expect(html).not.toContain('donor 0');
 expect(html).toContain('Missing observations are not PASS');
});
test('initiation observation is distinguished from full-CDS MFE',async({page})=>{
 await page.evaluate(()=>renderMfeWarning({metrics:{mfe_status:'not_computed',mfe_status_reason:'not_requested'},scientific_review:{initiation_mfe:{status:'AVAILABLE'}}}));
 await expect(page.locator('#mfeWarningBanner')).toHaveText(/Full-CDS RNA\/MFE.*separate initiation-window observation/);
});
test('ZIP dossier uses recorded version and checks rather than invented passes',async({page})=>{
 const dossier=await page.evaluate(()=>evidenceDossierMarkdown({product_version:'3.5.5',policy_decision:'BLOCK',qc_decision_matrix:[{criterion:'type_iis',result:'FAIL'}]}, {optimized_sequence:'ATGGCTTGA',metrics:{gc_percent:0,cai:0}},'Synthetic','2026-10-12'));
 expect(dossier).toContain('Product Version:** 3.5.5');
 expect(dossier).toContain('Recorded Type IIS check: FAIL');
 expect(dossier).toContain('Policy decision: BLOCK');
 expect(dossier).toContain('Initiation-window MFE: Not recorded');
 expect(dossier).toContain('GC Content | 0%');
 expect(dossier).not.toContain('cryptographic integrity manifest');
 expect(dossier).not.toContain('Zero Missense');
});

test('missing structure remains unreviewed and never requests auto resolution',async({page})=>{
 let autoCalls=0;
 page.on('request',request=>{if(request.postData()?.includes('auto_resolve_structure')) autoCalls++;});
 await page.locator('input[name="researchDesignMode"][value="paired_3dp"]').check();
 await page.locator('#sequenceInput').fill('MASTWQHLEKDPVNRSGYF');
 await page.locator('#recommendStructurePositions').click();
 await expect(page.locator('#structureWebReviewed')).not.toBeChecked();
 await expect(page.locator('#structureWebPositions')).toHaveValue('');
 expect(autoCalls).toBe(0);
});
test('viewer does not borrow a structure from another result',async({page})=>{
 const payload=await page.evaluate(()=>{
  state.structureViewerData={source_label:'stale',cartoon_pdb:'ATOM',ca_trace:[{aa:'M'}]};
  return buildStructureViewerPayload({input_type:'protein'}, {optimized_sequence:'ATG',metrics:{}});
 });
 expect(payload.source_label).toBe('NOT_AVAILABLE');
 expect(payload.cartoonPdb).toBeNull();
});
