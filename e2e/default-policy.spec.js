const {test,expect}=require('@playwright/test');
test('public UI keeps the default policy and ignores saved custom settings', async({page})=>{
 let payload;
 await page.addInitScript(()=>localStorage.setItem('factorforge_active_sop_v2',JSON.stringify({profile_id:'saved_custom',workflow:{watermark:true}})));
 await page.route('**/api/optimize',r=>{
  if(r.request().method()==='POST'){payload=r.request().postDataJSON();return r.fulfill({status:400,json:{error:'Synthetic capture only'}});}
  return r.fulfill({json:{supported_objectives:['dp_v2_1_1'],capabilities:{direct_structure:{available:true}},validation_checks:[]}});
 });
 await page.goto('/');
 await expect(page.locator('#sopProfileName')).toHaveText('FactorForge default policy');
 await expect(page.locator('#profileManagement')).toBeHidden();
 await expect(page.locator('#sopPresetSelect')).toBeHidden();
 await expect(page.locator('#manualSopOverrides')).toBeHidden();
 await expect(page.locator('#sopFileUpload')).toBeDisabled();
 await page.locator('#sequenceInput').fill('MASTWQHLEKDPVNRSGYF');
 await page.locator('#optimizeBtn').click();
 await expect.poll(()=>payload).toBeTruthy();
 expect(payload.sop_profile.profile_id).not.toBe('saved_custom');
 expect(payload.sop_profile.workflow.watermark).toBe(false);
 expect(payload.sop_profile.rules['assembly.type_iis.bsai.v1']).toBe('HARD_FAIL');
 const stored=await page.evaluate(()=>JSON.parse(localStorage.getItem('factorforge_active_sop_v2')));
 expect(stored.profile_id).toBe('saved_custom');
 await page.locator('input[name="researchDesignMode"][value="paired_3dp"]').check();
 await expect(page.locator('#structureWebInputs')).toBeVisible();
 await expect(page.locator('#appliedPolicySummary')).toContainText('3DP selected');
});
