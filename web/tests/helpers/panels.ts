// Shared Playwright helpers for the side panels, which both start collapsed on every page load.

import type { Page } from "@playwright/test";

async function openViaToolbar(page: Page, buttonName: string): Promise<void> {
  const button = page.getByRole("button", { name: buttonName, exact: true });
  // The toolbar is built by main.ts, so this also waits for the app to have mounted.
  await button.waitFor();
  if ((await button.getAttribute("aria-pressed")) !== "true") await button.click();
}

/** Opens the metrics (right) panel through its toolbar button. */
export async function openMetricsPanel(page: Page): Promise<void> {
  await openViaToolbar(page, "Panel");
}

/** Opens the parameters (left) panel through its toolbar button. */
export async function openParamsPanel(page: Page): Promise<void> {
  await openViaToolbar(page, "Parameters");
}
