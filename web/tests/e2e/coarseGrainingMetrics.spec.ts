import { expect, test } from "@playwright/test";
import { openField } from "../helpers/paramsPanel";

test.beforeEach(async ({ page }) => {
  await page.addInitScript(() => {
    localStorage.removeItem("milldynamics.paramsPanel");
    localStorage.removeItem("milldynamics.paramsGroups");
    localStorage.removeItem("milldynamics.panel");
  });
});

test("collision rate and impact histogram are hidden with an explanatory banner while coarse-graining is active, and reappear once k = 1", async ({
  page,
}) => {
  await page.goto("/");
  await page.locator(".hud", { hasText: /rpm/ }).waitFor();

  // Default load: 1 m drum / 10 mm balls / J = 0.30 gives N_true = 2460, well above the app's
  // actual boot max_balls (150, the Realtime quality preset -- worker.ts's INITIAL_PRESET_ID --
  // not mill-core's own SimulationParams::default() of 600), so k > 1 out of the box.
  const warnBanner = page.locator(".metrics-warn");
  await expect(warnBanner).toBeVisible();
  await expect(warnBanner).toHaveText(/Coarse-graining is active/);
  await expect(warnBanner).toContainText("2460");

  const collisionRateRow = page.locator(".metric-row", { has: page.locator(".metric-label", { hasText: "Collision rate" }) });
  await expect(collisionRateRow).toBeHidden();

  const histogramWrap = page.locator(".impact-histogram-wrap");
  await expect(histogramWrap).toBeHidden();

  // Raise the ball diameter enough that N_true drops below max_balls = 150 (N_true ~= 121 at
  // 45 mm, per derived.ts's trueBallCount formula), turning coarse-graining off (k = 1).
  const ballDiameterInput = await openField(page, "Ball diameter (mm)");
  await ballDiameterInput.fill("45");
  await page.getByRole("button", { name: "Apply", exact: true }).click();

  // Apply triggers a real WASM reinit (media.ball_diameter_m is resetRequired), which can take a
  // few seconds -- poll rather than using a fixed short timeout.
  await expect
    .poll(() => page.locator(".params-status").textContent(), { timeout: 15_000 })
    .toBe("Applied");

  await expect(warnBanner).toBeHidden();
  await expect(collisionRateRow).toBeVisible();
  await expect(histogramWrap).toBeVisible();
});
