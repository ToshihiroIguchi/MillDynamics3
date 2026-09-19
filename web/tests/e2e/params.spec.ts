import { expect, test } from "@playwright/test";

test("derived N_sim readout in the parameters modal never exceeds the max-balls input", async ({ page }) => {
  await page.goto("/");

  // The Parameters modal reads `state.params`, which only arrives once the worker's "ready"
  // message lands -- opening it before that leaves every derived-value row blank (see
  // ui/paramsModal.ts's `refreshDerived` early-return guard). Wait for the HUD's rpm readout
  // (also gated on `state.params`) so the click below is guaranteed to see real params.
  await page.locator(".hud", { hasText: /rpm/ }).waitFor();

  await page.getByRole("button", { name: "Parameters", exact: true }).click();
  const dialog = page.locator("dialog.params-modal");
  await expect(dialog).toBeVisible();

  const maxBallsInput = page
    .locator(".params-row", { has: page.locator("span", { hasText: "Max balls (coarse-graining target)" }) })
    .locator("input");
  const maxBallsValue = Number(await maxBallsInput.inputValue());

  const nSimDd = page.locator("dt", { hasText: "Simulated balls (N_sim)" }).locator("xpath=following-sibling::dd[1]");
  const nSimValue = Number(await nSimDd.textContent());

  expect(nSimValue).toBeLessThanOrEqual(maxBallsValue);

  // Stronger check: push N_true up by shrinking the ball diameter, and re-verify the invariant
  // still holds after refreshDerived() reacts to the input event.
  const ballDiameterInput = page
    .locator(".params-row", { has: page.locator("span", { hasText: "Ball diameter (mm)" }) })
    .locator("input");
  await ballDiameterInput.fill("5");
  await ballDiameterInput.dispatchEvent("input");

  const nSimValueAfter = Number(await nSimDd.textContent());
  const maxBallsValueAfter = Number(await maxBallsInput.inputValue());
  expect(nSimValueAfter).toBeLessThanOrEqual(maxBallsValueAfter);
});
