import { expect, test } from "@playwright/test";
import { openField } from "../helpers/paramsPanel";

test.beforeEach(async ({ page }) => {
  // Guarantee a clean, default panel/group state regardless of what a prior run/session left in
  // localStorage (panel visible, Mill/Media/Derived groups open, Lifters/Slurry/Simulation
  // collapsed -- see ui/paramsPanel.ts's DEFAULT_GROUP_OPEN).
  await page.addInitScript(() => {
    localStorage.removeItem("milldynamics.paramsPanel");
    localStorage.removeItem("milldynamics.paramsGroups");
  });
});

test("derived N_sim readout in the parameters panel never exceeds the max-balls input", async ({ page }) => {
  await page.goto("/");

  // The panel reads `state.params`, which only arrives once the worker's "ready" message lands --
  // reading it before that leaves every derived-value row blank (see ui/paramsPanel.ts's
  // `refreshDerived` early-return guard). Wait for the HUD's rpm readout (also gated on
  // `state.params`) so the reads below are guaranteed to see real params.
  await page.locator(".hud", { hasText: /rpm/ }).waitFor();

  await expect(page.locator("#params-panel")).toBeVisible();

  // "Max balls" lives in the "Simulation" group, which is collapsed by default.
  const maxBallsInput = await openField(page, "Max balls (coarse-graining target)");
  const maxBallsValue = Number(await maxBallsInput.inputValue());

  const nSimDd = page.locator("dt", { hasText: "Simulated balls (N_sim)" }).locator("xpath=following-sibling::dd[1]");
  const nSimValue = Number(await nSimDd.textContent());

  expect(nSimValue).toBeLessThanOrEqual(maxBallsValue);

  // Stronger check: push N_true up by shrinking the ball diameter (in the "Media" group, open by
  // default), and re-verify the invariant still holds after refreshDerived() reacts to the input
  // event.
  const ballDiameterInput = await openField(page, "Ball diameter (mm)");
  await ballDiameterInput.fill("5");
  await ballDiameterInput.dispatchEvent("input");

  const nSimValueAfter = Number(await nSimDd.textContent());
  const maxBallsValueAfter = Number(await maxBallsInput.inputValue());
  expect(nSimValueAfter).toBeLessThanOrEqual(maxBallsValueAfter);
});

test("applying an off-step value round-trips through the WASM core", async ({ page }) => {
  await page.goto("/");

  await page.locator(".hud", { hasText: /rpm/ }).waitFor();

  // Media density (schema step: 100) is open by default. 7850 is off-step -- with the old
  // <dialog>-based form's native HTML validation this silently failed to submit; the panel's
  // JS-side validation (form.noValidate = true) has no step check, so it must go through.
  const densityInput = await openField(page, "Media density (kg/m3)");
  await densityInput.fill("7850");

  // media.density_kg_m3 is `resetRequired` (schema.ts) -- baked into the ball population's mass
  // at seed time -- so this Apply still goes through a real WASM reinit, same as before Apply
  // started hot-applying fields that don't need one.
  await page.getByRole("button", { name: "Apply", exact: true }).click();

  // Apply triggers a real WASM reinit (resets the simulation), which can take a few seconds --
  // poll rather than using a fixed short timeout.
  await expect
    .poll(() => page.locator(".params-status").textContent(), { timeout: 15_000 })
    .toBe("Applied");
  await expect(page.locator(".params-status")).not.toHaveClass(/is-error/);

  // Confirms the value round-tripped through the Rust core and back via setParams(), not just
  // that the DOM was left untouched.
  await expect(densityInput).toHaveValue("7850");
});

test("applying a hot-appliable (resetRequired: false) field shows Applied, not a stuck Edited status", async ({
  page,
}) => {
  await page.goto("/");

  await page.locator(".hud", { hasText: /rpm/ }).waitFor();

  // Slurry viscosity (schema.ts: resetRequired: false) is in the "Slurry" group, collapsed by
  // default. This field hot-applies via "setParams" (no "ready" reply from the worker), which is
  // exactly the path that used to leave the status indicator stuck on "Edited -- not applied"
  // even though the apply succeeded (main.ts never called paramsPanel.setParams() on that path).
  const viscosityInput = await openField(page, "Viscosity (Pa·s)");
  await viscosityInput.fill("12.5");

  const status = page.locator(".params-status");
  await expect(status).toHaveText("Edited — not applied");

  await page.getByRole("button", { name: "Apply", exact: true }).click();

  await expect(status).toHaveText("Applied");
  await expect(status).toHaveClass(/is-applied/);
  await expect(status).not.toHaveClass(/is-edited/);
});

test("an out-of-range value shows a visible error without freezing the simulation", async ({ page }) => {
  await page.goto("/");

  await page.locator(".hud", { hasText: /rpm/ }).waitFor();

  // Drum diameter (schema min: 0.03) is in the "Mill" group, open by default.
  const diameterInput = await openField(page, "Drum diameter (m)");
  await diameterInput.fill("0");

  await page.getByRole("button", { name: "Apply", exact: true }).click();

  const errorBanner = page.locator(".params-error");
  await expect(errorBanner).toBeVisible();
  await expect(errorBanner).not.toHaveText("");

  const status = page.locator(".params-status");
  await expect(status).toHaveText("Error");
  await expect(status).toHaveClass(/is-error/);

  // The simulation must not be frozen by the failed Apply -- the HUD's rpm/frame readout should
  // keep advancing.
  const hud = page.locator(".hud", { hasText: /rpm/ });
  const hudTextBefore = await hud.textContent();
  await page.waitForTimeout(1000);
  await expect
    .poll(() => hud.textContent())
    .not.toBe(hudTextBefore);
});
