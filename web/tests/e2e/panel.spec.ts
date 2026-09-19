import { expect, test } from "@playwright/test";

test("metrics panel shows live values and toggling it resizes the canvas", async ({ page }) => {
  // Guarantee a clean, default-visible panel state regardless of what a prior run/session left in
  // localStorage.
  await page.addInitScript(() => localStorage.removeItem("milldynamics.panel"));
  await page.goto("/");

  await expect(page).toHaveTitle("MillDynamics3");
  await expect(page.locator("#scene")).toBeVisible();

  // Let the sim run for a few seconds so metrics have real values.
  await page.waitForTimeout(4000);

  const mixingIndexValue = page
    .locator(".metric-row", { has: page.locator(".metric-label", { hasText: "Mixing index" }) })
    .locator(".metric-value");
  await expect(mixingIndexValue).not.toHaveText("-");

  const widthBefore = await page.locator("#scene").evaluate((el) => el.clientWidth);

  await page.getByRole("button", { name: "Panel", exact: true }).click();
  await expect(page.locator(".metrics-panel")).toBeHidden();

  // The canvas resize is driven by a ResizeObserver reacting to the layout change, so it may land
  // a frame or two after the panel's visibility toggles -- poll rather than reading immediately.
  await expect
    .poll(() => page.locator("#scene").evaluate((el) => el.clientWidth))
    .toBeGreaterThan(widthBefore);
});
