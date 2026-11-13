const { test, expect } = require('@playwright/test');

test('release modal shows current public changes and explicit boundaries', async ({ page }) => {
  await page.goto('/');
  await page.locator('#changelogBtn').click();
  const modal = page.locator('#changelogModal');
  await expect(modal).toBeVisible();
  await expect(modal.locator('h3').first()).toHaveText('v3.6.8');
  await expect(modal).toContainText('Historical entries do not override');
  await expect(modal).toContainText('Public 1DP, 3DP, and 4DP Web Workflows');
  await expect(modal).toContainText('GNU AGPL-3.0');
  await expect(modal).toContainText('eijex.lab@gmail.com');
  await expect(modal).toContainText('not a guarantee of infrastructure-wide zero retention');
  await expect(modal).not.toContainText(/Internal Milestone|Stage 009|Training Telemetry|leak-free|67\/67/);
  await expect(modal.getByRole('link', {name:'Full software change history (GitHub)'}))
    .toHaveAttribute('href', 'https://github.com/eijex/factorforge-cds/blob/main/CHANGELOG.md');
});

test('release modal can be closed and reopened without retaining old content', async ({ page }) => {
  await page.goto('/');
  await page.locator('#changelogBtn').click();
  await page.locator('#closeModal').click();
  await expect(page.locator('#changelogModal')).toBeHidden();
  await page.locator('#changelogBtn').click();
  await expect(page.locator('#changelogModal h3').first()).toHaveText('v3.6.8');
});
