import { expect, test } from "@playwright/test";
import { openField } from "../helpers/paramsPanel";

test.beforeEach(async ({ page }) => {
  await page.addInitScript(() => {
    localStorage.removeItem("milldynamics.paramsPanel");
    localStorage.removeItem("milldynamics.paramsGroups");
    localStorage.removeItem("milldynamics.panel");
    localStorage.removeItem("milldynamics.metricsGroups");
  });
});

test("collision rate and impact histogram are hidden with an explanatory banner while coarse-graining is active, and reappear once k = 1", async ({
  page,
}) => {
  await page.goto("/");
  await page.locator(".hud", { hasText: /rpm/ }).waitFor();

  // The Grinding group (banner, Collision rate row, histogram) is collapsed by default -- the
  // "k = ..." chip on its summary is the only thing visible without opening it.
  const grindingGroup = page.locator(".metrics-group", { has: page.locator("summary", { hasText: "Grinding" }) });
  const grindingChip = grindingGroup.locator(".metrics-summary-chip");

  // Default load: 1 m drum / 10 mm balls / J = 0.30 gives N_true = 2460, well above the app's
  // actual boot max_balls (150, the Realtime quality preset -- worker.ts's INITIAL_PRESET_ID --
  // not mill-core's own SimulationParams::default() of 600), so k > 1 out of the box. The exact k
  // isn't asserted here (it depends on N_true/max_balls, not pinned to a round number) -- only
  // that the chip's format matches, and that it agrees with the banner's own k once opened.
  await expect(grindingChip).toBeVisible();
  await expect(grindingChip).toHaveText(/^k = \d+\.\d{2}$/);
  const chipText = await grindingChip.textContent();

  // Same "open a collapsed group programmatically" idiom as tests/helpers/paramsPanel.ts's
  // `openField` -- avoids any ambiguity about where a real click would land inside the summary
  // (which also contains the chip span above).
  await grindingGroup.evaluate((el) => {
    (el as HTMLDetailsElement).open = true;
  });

  const warnBanner = page.locator(".metrics-warn");
  await expect(warnBanner).toBeVisible();
  await expect(warnBanner).toHaveText(/Coarse-graining is active/);
  await expect(warnBanner).toContainText("2460");
  await expect(warnBanner).toContainText(`(${chipText})`);

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

  await expect(grindingChip).toBeHidden();
  await expect(warnBanner).toBeHidden();
  await expect(collisionRateRow).toBeVisible();
  await expect(histogramWrap).toBeVisible();
});
