// @ts-check
const { test, expect } = require('@playwright/test');
const path = require('node:path');
const { pathToFileURL } = require('node:url');

test.use({ deviceScaleFactor: 2 });

test.describe('HPX Evidence Viewer Offline Browser Regression', () => {
  test('reviewed-index.html renders with networking disabled with zero external network requests and active Canvas interactions', async ({ page, context }, testInfo) => {
    const networkRequests = [];
    const pageErrors = [];
    page.on('pageerror', error => pageErrors.push(error.message));
    page.on('request', (req) => {
      networkRequests.push(req.url());
    });

    // 1. Force strict offline mode
    await context.setOffline(true);

    // 2. Load reviewed viewer via file protocol
    const viewerPath = path.resolve(__dirname, '../examples/hpx_evidence_viewer/reviewed-index.html');
    const fileUrl = pathToFileURL(viewerPath).href;

    await page.goto(fileUrl, { waitUntil: 'load' });

    // Assert zero remote/external network requests
    const externalRequests = networkRequests.filter((url) => !url.startsWith('file://') && !url.startsWith('data:'));
    expect(externalRequests).toEqual([]);

    // 3. Verify 3D canvas exists and has dimensions
    const cv3D = page.locator('#canvas3d');
    await expect(cv3D).toBeVisible();
    const box3D = await cv3D.boundingBox();
    expect(box3D).not.toBeNull();
    expect(box3D.width).toBeGreaterThan(200);
    expect(box3D.height).toBeGreaterThan(200);

    // 4. Verify track canvases exist and have dimensions
    for (const id of ['#cvCoord', '#cvDomains', '#cvPlddt', '#cvDesign', '#cvAdapt']) {
      const loc = page.locator(id);
      await expect(loc).toBeVisible();
      const b = await loc.boundingBox();
      expect(b.width).toBeGreaterThan(300);
      expect(b.height).toBeGreaterThan(20);
    }

    // 5. Verify initial provenance & inspector state
    await expect(page.locator('#inspWi')).toHaveText('NOT_AVAILABLE');
    await expect(page.locator('#inspRes')).toContainText('240');

    // 6. Test interaction: Click on coordinates canvas (e.g. residue 50 area ~ 10% across)
    const cvCoord = page.locator('#cvCoord');
    const coordBox = await cvCoord.boundingBox();
    expect(coordBox).not.toBeNull();

    // Click near residue 50
    await page.mouse.click(coordBox.x + coordBox.width * 0.11, coordBox.y + coordBox.height / 2);
    await page.waitForTimeout(100);

    // Selected residue should be updated in inspector
    const resText = await page.locator('#inspRes').textContent();
    expect(resText).toMatch(/#51\s+[A-Za-z]{3}/i);

    // 7. Test Construct switcher (Precursor vs Mature)
    await page.locator('#btnMature').click();
    await expect(page.locator('#inspPos')).toContainText('28 (Construct)');

    // 8. Test responsive resize stability
    const originalHeight = coordBox.height;
    await page.setViewportSize({ width: 900, height: 700 });
    await page.waitForTimeout(100);
    const resizedBox = await cvCoord.boundingBox();
    expect(resizedBox.height).toBe(originalHeight);

    // 9. Capture screenshot to artifact directory as visual evidence
    const screenshotPath = testInfo.outputPath('hpx-reviewed-viewer.png');
    await page.screenshot({ path: screenshotPath, fullPage: true });
    expect(pageErrors).toEqual([]);
  });

  test('canonical index.html renders with networking disabled with zero external network requests', async ({ page, context }, testInfo) => {
    const networkRequests = [];
    const pageErrors = [];
    page.on('pageerror', error => pageErrors.push(error.message));
    page.on('request', (req) => {
      networkRequests.push(req.url());
    });

    await context.setOffline(true);
    const viewerPath = path.resolve(__dirname, '../examples/hpx_evidence_viewer/index.html');
    const fileUrl = pathToFileURL(viewerPath).href;

    await page.goto(fileUrl, { waitUntil: 'load' });

    const externalRequests = networkRequests.filter((url) => !url.startsWith('file://') && !url.startsWith('data:'));
    expect(externalRequests).toEqual([]);

    // Check canvases in index.html
    const cvCoord = page.locator('#cvCoord');
    await expect(cvCoord).toBeVisible();
    expect(pageErrors).toEqual([]);
  });
});
