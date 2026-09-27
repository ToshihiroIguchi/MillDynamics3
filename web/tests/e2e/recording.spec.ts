import { expect, test } from "@playwright/test";

test("Record toggles a recording and downloads a WebM video on stop", async ({ page }) => {
  await page.addInitScript(() => localStorage.removeItem("milldynamics.panel"));
  await page.goto("/");

  const recordBtn = page.getByRole("button", { name: "Record", exact: true });
  await expect(recordBtn).toBeEnabled();
  await recordBtn.click();

  await expect(page.getByRole("button", { name: "Stop recording" })).toBeVisible();

  // A couple of animation frames' worth of capture is enough for a smoke test -- this isn't
  // checking capture quality, just that start -> stop -> download works end to end.
  await page.waitForTimeout(500);

  const downloadPromise = page.waitForEvent("download");
  await page.getByRole("button", { name: "Stop recording" }).click();
  const download = await downloadPromise;

  expect(download.suggestedFilename()).toMatch(/^milldynamics-recording-.*\.webm$/);
});
