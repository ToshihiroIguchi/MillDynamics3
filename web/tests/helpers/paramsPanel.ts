// Shared Playwright helpers for driving the params panel (ui/paramsPanel.ts) from e2e specs.
// Lives outside tests/e2e (playwright.config.ts's testDir) so Playwright doesn't try to run this
// file itself as a spec.

import type { Page } from "@playwright/test";

/** Locates a field's `<input>` by its row's label text, forcing the row's group `<details>` open
 * first (fields in collapsed groups aren't fillable/readable via Playwright). */
export async function openField(page: Page, labelText: string) {
  const row = page.locator(".params-row", { has: page.locator("span", { hasText: labelText }) });
  // Every ancestor <details> (the group and, for rarely-used fields, its Advanced sub-section).
  const ancestors = row.locator("xpath=ancestor::details");
  const count = await ancestors.count();
  for (let i = 0; i < count; i++) {
    await ancestors.nth(i).evaluate((el) => {
      (el as HTMLDetailsElement).open = true;
    });
  }
  return row.locator("input");
}
