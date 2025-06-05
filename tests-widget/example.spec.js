// @ts-check
import { test, expect } from '@playwright/test';

test('Renders Vitessce widget containing a scatterplot view (Jupyter)', async ({ page }) => {
  test.setTimeout(60_000);
  await page.goto('http://localhost:3000/widget_from_dict.html');

  // Expect a title "to contain" a substring.
  await expect(page).toHaveTitle('widget_from_dict');

  await expect(page.getByText('Scatterplot (UMAP)')).toBeVisible();
  await expect(page.getByText('523 cells')).toHaveCount(3);
});

test('Renders Vitessce widget containing a scatterplot view (Marimo)', async ({ page }) => {
  test.setTimeout(60_000);
  await page.goto('http://localhost:3000/marimo.html');

  // Expect a title "to contain" a substring.
  await expect(page).toHaveTitle('marimo');

  await expect(page.getByText('Scatterplot (UMAP)')).toBeVisible();
  await expect(page.getByText('523 cells')).toHaveCount(3);
});